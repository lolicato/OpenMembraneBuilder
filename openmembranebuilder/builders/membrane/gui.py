"""The widgets of the membrane builder: what you see and fill in on its pages.

- **Membrane page** (`render_membrane`): force field, box, solvation, number of systems,
  output name in the sidebar; template choice, lipid table and composition rows in the main
  area.
- **Protein insertion page** (`render_protein_insertion`): the box again, then for every
  protein its position, height, rotation and copies in the sidebar; the schematic system view
  with warnings in the main area.
- **Preview page** (`render_preview`): a summary of all settings and the Build button.

The values are kept in Streamlit's session memory under fixed names (keys), so they survive
the redraw after every click and going back and forth between pages.
"""

import copy
import html
import os
import re
from pathlib import Path

import pandas as pd
import streamlit as st

# The local duplicate of GlobalGui was removed; use the single shared one.
# (src/core/gui.py must not call st.set_page_config any more - it is called in app.py.)
from openmembranebuilder.core.gui import GlobalGui


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _ns(name: str) -> str:
    """A short prefix for the widget names of one protein, made from its file name.

    Args:
        name: The protein's file name, e.g. `"1UBQ.pdb"`.

    Returns:
        The prefix, e.g. `"1UBQ_pdb"`.
    """
    return name.replace(".", "_").replace(" ", "_")


def _default(d: dict, name: str, fallback):
    """A saved value, or a fallback when it is missing or empty.

    Args:
        d: The saved values.
        name: The value's name.
        fallback: What to use instead.

    Returns:
        `d[name]`, or `fallback` when it is missing, `None` or empty.
    """
    val = d.get(name)
    return fallback if val is None or val == "" else val


def _num_input(label: str, key: str, lo: float, hi: float, step: float, default: float, container=None):
    """A number field that always shows a value inside its allowed range.

    The stored value is set before the field is drawn and kept within `[lo, hi]`, so
    making the box smaller can never leave a position outside the box.

    Args:
        label: Text next to the field.
        key: Name under which the value is remembered.
        lo: Smallest allowed value.
        hi: Largest allowed value.
        step: Change per click on the +/− buttons.
        default: Value used the first time.
        container: Where to draw the field; the sidebar by default.

    Returns:
        The value shown in the field.
    """
    if key not in st.session_state:
        st.session_state[key] = default
    st.session_state[key] = min(max(float(st.session_state[key]), float(lo)), float(hi))
    return (container or st.sidebar).number_input(
        label, min_value=float(lo), max_value=float(hi), step=float(step), key=key
    )


def _reset_entry_widgets() -> None:
    """Forget the values of the lipid-row widgets, so they are filled again from the rows.

    Needed after a row is deleted or a template is loaded: the rows move, and the widgets
    would otherwise keep showing the values of the old rows.
    """
    for k in list(st.session_state.keys()):
        if re.match(r"^lipid_\d+$", k) or re.match(r"^entry_\d+_\d+$", k):
            del st.session_state[k]


def _detect_template_abs_mode(template_path: str) -> bool | None:
    """Guess whether a template gives numbers of lipids or fractions.

    The first lipid row decides: whole numbers mean numbers of lipids, anything else
    fractions.

    Args:
        template_path: The template file.

    Returns:
        `True` for numbers, `False` for fractions, `None` when the file cannot be read or
        has no lipid row.
    """
    try:
        with open(template_path, "r", encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None

    for raw in lines:
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        cols = [part.strip() for part in [p.strip() for p in line.split(",")]]
        if len(cols) < 3:
            continue
        try:
            upper = float(cols[1])
            lower = float(cols[2])
        except ValueError:
            continue

        is_absolute = abs(upper - round(upper)) < 1e-9 and abs(lower - round(lower)) < 1e-9
        return is_absolute
    return None


# ─── Martini version consistency ──────────────────────────────────────────────

def _ff_major(ff: str) -> str | None:
    """The Martini generation of a force-field name.

    Args:
        ff: A name such as `"martini_v3"`, `"martini3001"` or `"martini22"`.

    Returns:
        `"2"`, `"3"`, or `None` when it cannot be told.
    """
    ff = (ff or "").lower()
    if "v2" in ff or "martini2" in ff:
        return "2"
    if "v3" in ff or "martini3" in ff:
        return "3"
    return None


def ff_mismatch_message(protein_files: dict, selected_ff: str) -> str | None:
    """A warning when proteins and membrane belong to different Martini generations.

    Martini 2 and Martini 3 models must not be mixed: the simulation may fail or give
    wrong results.

    Args:
        protein_files: The proteins, each with its `force_field`.
        selected_ff: The force field chosen for the membrane.

    Returns:
        The warning text, or `None` when everything matches or is unknown.
    """
    versions = {
        name: _ff_major(p.get("force_field", ""))
        for name, p in (protein_files or {}).items()
        if _ff_major(p.get("force_field", ""))
    }
    if not versions:
        return None

    parts = []
    if len(set(versions.values())) > 1:
        listing = ", ".join(f"{n}: Martini {v}" for n, v in versions.items())
        parts.append(f"The proteins use different Martini versions ({listing}).")

    sel = _ff_major(selected_ff)
    if sel:
        wrong = [f"{n} (Martini {v})" for n, v in versions.items() if v != sel]
        if wrong:
            parts.append(
                f"The membrane force field is {selected_ff} (Martini {sel}), "
                f"but these proteins were made for another version: {', '.join(wrong)}."
            )

    if not parts:
        return None
    return " ".join(parts) + (
        " Martini 2 and Martini 3 models are not compatible; mixing them can make the "
        "simulation fail or give physically wrong results."
    )


# ─── Membrane page widgets ────────────────────────────────────────────────────

class MembraneGui:
    """The widgets of the membrane page."""

    @staticmethod
    def render_force_field_selector(force_fields: list[str]) -> str:
        """The force-field choice in the sidebar (Martini 3 by default).

        Args:
            force_fields: Names of the available force fields.

        Returns:
            The chosen force field.
        """
        if "selected_ff" not in st.session_state or st.session_state["selected_ff"] not in force_fields:
            st.session_state["selected_ff"] = "martini_v3" if "martini_v3" in force_fields else force_fields[0]
        return st.sidebar.selectbox("Force Field", force_fields, key="selected_ff")

    @staticmethod
    def render_box_params() -> tuple:
        """Box type and size in the sidebar (also used on the protein insertion page).

        Only the lengths that the box type needs are shown: x, y and z for a rectangular box,
        x and z for a hexagonal box, x for a skewed hexagonal box.

        Returns:
            `(box type, x, y, z)`, lengths in nm.
        """
        st.sidebar.subheader("Box Parameters")
        st.session_state.setdefault("box_type", "rectangular")
        st.session_state.setdefault("boxx", 10.0)
        st.session_state.setdefault("boxy", 10.0)
        st.session_state.setdefault("boxz", 10.0)

        selected_box_type = st.sidebar.selectbox(
            "Box Type",
            ["rectangular", "hexagonal", "skewed_hexagonal"],
            key="box_type",
            help="rectangular: X,Y,Z | hexagonal: X,Z | skewed_hexagonal: X only",
        )

        col1, col2, col3 = st.sidebar.columns(3)
        boxx = col1.number_input("Box X (nm)", min_value=5.0, max_value=50.0, key="boxx")
        boxy = (
            col2.number_input("Box Y (nm)", min_value=5.0, max_value=50.0, key="boxy")
            if selected_box_type == "rectangular"
            else st.session_state.get("boxy", 10.0)
        )
        boxz = (
            col3.number_input("Box Z (nm)", min_value=5.0, max_value=50.0, key="boxz")
            if selected_box_type != "skewed_hexagonal"
            else st.session_state.get("boxz", 10.0)
        )
        return selected_box_type, boxx, boxy, boxz

    @staticmethod
    def render_solvation_input(selected_ff: str) -> str:
        """Ions and salt concentration in the sidebar.

        With Martini 2.2 a share of the water can be replaced by antifreeze water (WF).

        Args:
            selected_ff: The chosen force field.

        Returns:
            The solvation in COBY's notation, e.g. `"solv:W pos:NA neg:CL salt_molarity:0.15"`.
        """
        st.sidebar.subheader("Solvation")
        st.session_state.setdefault("pos_ion", "NA")
        st.session_state.setdefault("neg_ion", "CL")
        st.session_state.setdefault("salt_molarity", 0.15)
        st.session_state.setdefault("replace_w_percentage", 0.0)

        pos_ion = st.sidebar.selectbox("Positive Ion", ["NA", "CA"], key="pos_ion")
        neg_ion = st.sidebar.selectbox(
            "Negative Ion",
            ["CL", "BR", "IOD", "ACE", "BF4", "PF6", "SCN", "CLO4", "NO3"],
            key="neg_ion",
        )
        salt_molarity = st.sidebar.number_input(
            "Salt Molarity (M)", min_value=0.0, max_value=2.0, key="salt_molarity"
        )

        if selected_ff == "martini_v2.2":
            percentage_w_replace = st.sidebar.number_input(
                "Percentage of unfreezable water",
                min_value=0.0,
                max_value=100.0,
                key="replace_w_percentage",
                help="Percentage of water to replace with unfreezable water (WF)",
            )
            percentage_w = 100 - percentage_w_replace
            if percentage_w_replace != 0.0:
                return (
                    f"solv:W{percentage_w} solv:WF{percentage_w_replace} "
                    f"params:IMPORTED pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt_molarity}"
                )

        return f"solv:W pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt_molarity}"

    @staticmethod
    def render_n_systems_input() -> int:
        """The number of independent systems (replicas) to build, 1 to 99.

        Returns:
            The chosen number.
        """
        st.session_state.setdefault("n_systems", 1)
        return st.sidebar.number_input(
            "Systems",
            min_value=1,
            max_value=99,
            key="n_systems",
            help="Number of independent systems to create",
        )

    @staticmethod
    def render_output_name_input() -> str:
        """An optional name that is added to the output folder.

        Returns:
            The name typed, or an empty text.
        """
        st.session_state.setdefault("output_name", "")
        return st.sidebar.text_input("Output folder name", key="output_name")

    @staticmethod
    def render_predefined_template_selector(available_lipids: list, global_info):
        """A list of the predefined membrane templates (older single-list version of the choice).

        Args:
            available_lipids: The lipids of the chosen force field.
            global_info: The standard folder locations.

        Returns:
            The path of the chosen template, or `None` for "No template".
        """
        template_dir = Path(global_info.resources_path) / "membrane_templates"
        template_files = sorted(
            [p for p in template_dir.iterdir() if p.is_file() and p.suffix.lower() in {".csv", ".ob"}],
            key=lambda p: p.name.lower(),
        ) if template_dir.exists() else []

        if not template_files:
            return None

        template_source = st.session_state.setdefault(
            "template_source",
            {"kind": "manual", "name": "", "path": "", "loaded_for": "", "entries": [], "templates": {}},
        )

        options = ["No template"] + [p.stem for p in template_files]
        current_name = template_source.get("name", "")
        default_index = 0 if current_name not in options else options.index(current_name)
        selected_name = st.selectbox(
            "Predefined membrane template",
            options,
            index=default_index,
            key="builtin_template_selector",
        )

        if selected_name == "No template":
            if template_source.get("kind") == "upload" and template_source.get("path"):
                return template_source.get("path")
            if template_source.get("kind") != "manual":
                _reset_entry_widgets()
                backup = template_source.get("manual_backup")
                restored = copy.deepcopy(backup) if isinstance(backup, list) else []
                template_source.update({
                    "kind": "manual",
                    "name": "",
                    "path": "",
                    "loaded_for": "",
                    "entries": restored,
                    "templates": {},
                })
                st.session_state["entries"] = restored
            return None

        selected_path = str(next(p for p in template_files if p.stem == selected_name))
        if (
            template_source.get("path") != selected_path
            or template_source.get("kind") != "builtin"
            or template_source.get("loaded_for") != selected_path
        ):
            if template_source.get("kind") == "manual":
                current_entries = st.session_state.get("entries")
                if isinstance(current_entries, list):
                    template_source["manual_backup"] = copy.deepcopy(current_entries)
            detected_abs = _detect_template_abs_mode(selected_path)
            if detected_abs is not None:
                st.session_state["abs_lip_vals"] = detected_abs
                template_source["abs_lip_vals"] = detected_abs
            _reset_entry_widgets()
            template_source.update({
                "kind": "builtin",
                "name": selected_name,
                "path": selected_path,
                "loaded_for": "",
                "entries": [],
                "templates": {},
            })

        return selected_path

    @staticmethod
    def render_template_uploader(
        available_lipids: list,
        global_info,
        temp_folder,
        template_source: dict | None = None,
        entries: list | None = None,
        abs_lip_vals: bool | None = None,
    ) -> bool:
        """The template box: no template, a predefined template, or your own file.

        A template fills the lipid table. The rows typed by hand are kept aside and come back
        when "No template" is chosen again. An example template can be downloaded.

        Args:
            available_lipids: The lipids of the chosen force field (to check the template).
            global_info: The standard folder locations.
            temp_folder: Where an uploaded template is stored.
            template_source: The current template state.
            entries: The rows shown at the moment.
            abs_lip_vals: `True` when the rows give numbers of lipids.

        Returns:
            `True` when a template is in use.
        """
        os.makedirs(temp_folder, exist_ok=True)

        template_source = template_source or st.session_state.setdefault(
            "template_source",
            {"kind": "manual", "name": "", "path": "", "loaded_for": "", "entries": [], "templates": {}, "manual_backup": [], "abs_lip_vals": False},
        )
        template_active = template_source.get("kind") in {"builtin", "upload"}
        template_path = template_source.get("path", "")

        with st.expander("Do you have a membrane template or want to use a predefined one?", expanded=True):
            st.markdown(
                "Choose whether to use a predefined membrane template or upload your own "
                "(`.ob` or `.csv`) to automatically load membrane compositions."
            )

            current_kind = template_source.get("kind", "manual")
            source_choice = st.radio(
                "Template source",
                ["No template", "Use predefined template", "Upload template"],
                index=0 if current_kind == "manual" else 1 if current_kind == "builtin" else 2,
                horizontal=True,
                key="template_source_choice",
            )

            def apply_template_state(kind: str, name: str, path: str, entries: list | None, abs_mode: bool | None):
                """Switch the lipid table to a template's rows, or back to the rows typed by hand.
                """
                current_path = template_source.get("path", "")
                if current_kind in {"builtin", "upload"} and current_path and current_path != path:
                    if current_kind == "upload" and os.path.exists(current_path):
                        try:
                            os.remove(current_path)
                        except OSError:
                            pass
                backup = template_source.get("manual_backup")
                if current_kind == "manual" and isinstance(st.session_state.get("entries"), list):
                    template_source["manual_backup"] = copy.deepcopy(st.session_state["entries"])
                elif isinstance(backup, list):
                    template_source["manual_backup"] = copy.deepcopy(backup)
                if kind == "manual":
                    restored = copy.deepcopy(template_source.get("manual_backup", [])) if isinstance(template_source.get("manual_backup"), list) else []
                    template_source.update({
                        "kind": "manual",
                        "name": "",
                        "path": "",
                        "loaded_for": "",
                        "entries": restored,
                        "templates": {},
                        "abs_lip_vals": bool(abs_mode) if abs_mode is not None else bool(template_source.get("abs_lip_vals", False)),
                    })
                    st.session_state["entries"] = list(restored)
                    st.session_state["abs_lip_vals"] = bool(template_source.get("abs_lip_vals", False))
                    _reset_entry_widgets()
                    return

                parsed_entries = list(entries or [])
                abs_value = bool(abs_mode) if abs_mode is not None else bool(template_source.get("abs_lip_vals", False))
                template_source.update({
                    "kind": kind,
                    "name": name,
                    "path": path,
                    "loaded_for": path,
                    "entries": parsed_entries,
                    "templates": {},
                    "abs_lip_vals": abs_value,
                })
                st.session_state["entries"] = list(parsed_entries)
                st.session_state["abs_lip_vals"] = abs_value
                _reset_entry_widgets()

            if source_choice == "No template":
                if current_kind != "manual":
                    apply_template_state("manual", "", "", template_source.get("manual_backup", []), None)
                return False

            if source_choice == "Use predefined template":
                template_dir = Path(global_info.resources_path) / "membrane_templates"
                predefined_templates = sorted(
                    [p for p in template_dir.iterdir() if p.is_file() and p.suffix.lower() in {".csv", ".ob"}],
                    key=lambda p: p.name.lower(),
                ) if template_dir.exists() else []

                if not predefined_templates:
                    st.info("No predefined membrane templates are available in the resources folder.")
                    return False

                options = [p.stem for p in predefined_templates]
                current_name = template_source.get("name", "")
                default_index = 0 if current_name not in options else options.index(current_name)
                selected_name = st.selectbox(
                    "Predefined membrane template",
                    options,
                    index=default_index,
                    key="builtin_template_selector",
                )

                selected_path = str(next(p for p in predefined_templates if p.stem == selected_name))
                if template_source.get("kind") == "builtin" and template_source.get("path") == selected_path:
                    return True
                detected_abs = _detect_template_abs_mode(selected_path)
                parser = __import__("openmembranebuilder.builders.membrane.parser", fromlist=["MembraneParser"]).MembraneParser()
                template_data = parser.load_membrane_template_from_csv(selected_path, available_lipids, gui=False) or {}
                loaded_entries = list(next(iter(template_data.values()))) if template_data else []
                if not loaded_entries:
                    return False
                apply_template_state("builtin", selected_name, selected_path, loaded_entries, detected_abs)
                return True

            if source_choice == "Upload template":
                example_path = os.path.join(global_info.membrane_template_example_path, "example_template.ob")
                if os.path.exists(example_path):
                    with open(example_path, "rb") as f:
                        st.download_button(
                            "Download example template",
                            data=f.read(),
                            file_name="example_template.ob",
                        )

                uploaded_template = st.file_uploader("Template file", type=["csv", "ob"], label_visibility="collapsed")
                if uploaded_template:
                    template_path = os.path.join(temp_folder, uploaded_template.name)
                    if template_source.get("kind") == "upload" and template_source.get("path") == template_path:
                        return True
                    with open(template_path, "wb") as f:
                        f.write(uploaded_template.getvalue())
                    detected_abs = _detect_template_abs_mode(template_path)
                    parser = __import__("openmembranebuilder.builders.membrane.parser", fromlist=["MembraneParser"]).MembraneParser()
                    template_data = parser.load_membrane_template_from_csv(template_path, available_lipids, gui=False) or {}
                    loaded_entries = list(next(iter(template_data.values()))) if template_data else []
                    if not loaded_entries:
                        return False
                    apply_template_state("upload", uploaded_template.name, template_path, loaded_entries, detected_abs)
                    return False

        return False

    @staticmethod
    def render_lipid_mapping(lipid_map: dict):
        """The searchable table of all lipids of the force field, shown on request.

        It lists head group, linker, tails and bead types of every lipid and can be filtered
        and downloaded as CSV.

        Args:
            lipid_map: The lipid descriptions of the chosen force field.
        """
        if "show_lipid_mapping" not in st.session_state:
            st.session_state.show_lipid_mapping = False

        label = "Hide Lipid Mapping" if st.session_state.show_lipid_mapping else "Show Lipid Mapping"
        if st.button(label, key="toggle_lipid_mapping"):
            st.session_state.show_lipid_mapping = not st.session_state.show_lipid_mapping
            st.rerun()

        if not st.session_state.show_lipid_mapping or not lipid_map:
            return

        df = pd.DataFrame.from_dict(lipid_map, orient="index").reset_index(drop=True)
        st.subheader("Lipid Mapping")

        filter_cols = [c for c in df.columns if c != "resname"]
        cols = st.columns(len(filter_cols) + 1)
        active_filters = {}

        st.session_state.setdefault("lipid_map_resname_search", "")
        cols[0].text_input("Resname", placeholder="e.g. POPC", key="lipid_map_resname_search")

        for i, col in enumerate(filter_cols):
            options = ["All"] + sorted(df[col].dropna().astype(str).unique().tolist())
            key = f"lipid_filter_{col}"
            if key not in st.session_state or st.session_state[key] not in options:
                st.session_state[key] = "All"
            active_filters[col] = cols[i + 1].selectbox(col.capitalize(), options, key=key)

        filtered = df.copy()
        resname_search = st.session_state.get("lipid_map_resname_search", "")
        if resname_search:
            filtered = filtered[
                filtered["resname"].astype(str).str.contains(resname_search, case=False, na=False, regex=False)
            ]
        for col, val in active_filters.items():
            if val != "All":
                filtered = filtered[filtered[col].astype(str) == val]

        st.caption(f"Showing {len(filtered)} of {len(df)} lipids")
        st.dataframe(
            filtered,
            use_container_width=True,
            hide_index=True,
            column_config={
                "resname": st.column_config.TextColumn("Resname"),
                "moltype": st.column_config.TextColumn("Type"),
                "head": st.column_config.TextColumn("Headgroup"),
                "linker": st.column_config.TextColumn("Linker"),
                "tail1": st.column_config.TextColumn("Tail 1"),
                "tail2": st.column_config.TextColumn("Tail 2"),
                "tail1_beads": st.column_config.TextColumn("Tail 1 Beads"),
                "tail2_beads": st.column_config.TextColumn("Tail 2 Beads"),
            },
        )
        st.download_button(
            "Download as CSV",
            data=filtered.to_csv(index=False).encode("utf-8"),
            file_name="lipid_mapping.csv",
            mime="text/csv",
            key="lipid_map_download",
        )

    @staticmethod
    def render_membrane_entries(
        available_lipids: list,
        entries: list | None = None,
        abs_lip_vals: bool | None = None,
        template_source: dict | None = None,
    ) -> list:
        """The lipid table: one row per lipid, with its share in each leaflet.

        Columns: lipid, upper leaflet (RU/NU), lower leaflet (RL/NL), area per lipid upper
        (AU) and lower (AL), in nm². With "Give absolute lipid values" the leaflet columns
        are numbers of lipids (NU/NL), otherwise fractions (RU/RL).

        Args:
            available_lipids: The lipids of the chosen force field.
            entries: The rows to show.
            abs_lip_vals: `True` for numbers of lipids.
            template_source: The current template state.

        Returns:
            The rows as edited, `[lipid, upper, lower, AU, AL]` each.
        """
        default_abs = bool(st.session_state.get("abs_lip_vals", False if abs_lip_vals is None else abs_lip_vals))
        if abs_lip_vals is not None:
            st.session_state.setdefault("abs_lip_vals", default_abs)
        abs_lip_vals = bool(
            st.checkbox(
                "Give absolute lipid values",
                key="abs_lip_vals",
                value=default_abs,
            )
        )

        template_source = template_source or st.session_state.get("template_source", {
            "kind": "manual", "name": "", "path": "", "loaded_for": "", "entries": [], "templates": {}
        })
        if template_source.get("kind") in {"builtin", "upload"}:
            st.info("Lipid entries are loaded from the template.")
            current_path = template_source.get("path", "")
            current_entries = st.session_state.get("entries")
            if isinstance(current_entries, list) and current_entries != template_source.get("entries", []):
                template_source["entries"] = list(current_entries)
                st.session_state["template_source"] = template_source
            if current_path and (template_source.get("loaded_for") != current_path or not template_source.get("entries")):
                incoming = list(template_source.get("entries") or [])
                template_source["loaded_for"] = current_path
                template_source["entries"] = incoming
                st.session_state["entries"] = incoming
                _reset_entry_widgets()
                st.session_state["template_source"] = template_source

        st.markdown("**Membrane Entries**")

        if entries is None:
            entries = st.session_state.get("entries")
        if not isinstance(entries, list) or not entries:
            entries = [["POPC", 1.0, 1.0, 0.6, 0.6]]
        if st.session_state.get("entries") != list(entries):
            st.session_state["entries"] = list(entries)
            _reset_entry_widgets()

        def add_entry():
            """Add a new row (first lipid of the list, 1.0 / 1.0, area per lipid 0.6 nm²)."""
            st.session_state.entries.append(
                [available_lipids[0] if available_lipids else "POPC", 1.0, 1.0, 0.6, 0.6]
            )

        def delete_entry(idx):
            """Remove one row from the lipid table."""
            entries_now = st.session_state.get("entries")
            if not isinstance(entries_now, list):
                return
            updated = [entry for i, entry in enumerate(entries_now) if i != idx]
            st.session_state["entries"] = updated
            if isinstance(template_source, dict):
                template_source["entries"] = updated
                if template_source.get("kind") in {"builtin", "upload"}:
                    template_source["loaded_for"] = template_source.get("path", "")
            _reset_entry_widgets()
            st.session_state["_membrane_dirty"] = True

        headers = ["Lipid", "NU" if abs_lip_vals else "RU", "NL" if abs_lip_vals else "RL", "AU", "AL"]
        cols = st.columns([2, 1, 1, 1, 1, 1])
        for i, h in enumerate(headers):
            cols[i].write(h)

        for idx, entry in enumerate(st.session_state.entries):
            cols = st.columns([2, 1, 1, 1, 1, 1])

            lipid_key = f"lipid_{idx}"
            if lipid_key not in st.session_state or st.session_state[lipid_key] not in available_lipids:
                st.session_state[lipid_key] = entry[0] if entry[0] in available_lipids else available_lipids[0]
            with cols[0]:
                entry[0] = st.selectbox(f"Lipid {idx}", available_lipids, key=lipid_key)

            limits = (0, 10000, 1.0) if abs_lip_vals else (0.0, 1.0, 0.01)
            for i in range(1, 5):
                with cols[i]:
                    mn, mx, step = limits
                    state_key = f"entry_{idx}_{i}"
                    if state_key not in st.session_state:
                        st.session_state[state_key] = float(entry[i])
                    # switching absolute <-> relative changes the allowed range
                    st.session_state[state_key] = min(max(float(st.session_state[state_key]), float(mn)), float(mx))
                    entry[i] = st.number_input(
                        headers[i],
                        min_value=float(mn),
                        max_value=float(mx),
                        step=step,
                        key=state_key,
                    )
            with cols[5]:
                if st.button("🗑", key=f"del_{idx}"):
                    delete_entry(idx)
                    st.rerun()

        if st.button("Add lipid"):
            add_entry()
            st.rerun()
        return st.session_state.entries


# ─── Protein page widgets ─────────────────────────────────────────────────────

class ProteinGui:
    """The widgets of the protein insertion page."""

    @staticmethod
    def render_protein_list() -> dict[str, dict]:
        """The list of prepared proteins in the sidebar (read only).

        The proteins are changed in the protein-preparation step, not here.

        Returns:
            The proteins: `{file name: {"pdb_path", "itp_path", "force_field"}}`.
        """
        protein_files = st.session_state.get("protein_files", {})

        st.sidebar.subheader("Proteins")
        if not protein_files:
            st.sidebar.warning("No proteins available. Go back to the martinize builder and add a protein.")
        for pdb_name, paths in protein_files.items():
            st.sidebar.markdown(f"**✅ {pdb_name}**")
            st.sidebar.caption(f"ITP: {os.path.basename(paths['itp_path'])}")
        if protein_files:
            st.sidebar.caption("To change the proteins, go back to the martinize builder.")

        return protein_files

    @staticmethod
    def render_insertion_params_for_protein(
        protein_name: str,
        boxx: float,
        boxy: float,
        boxz: float,
        n_systems: int,
        protein_only: bool,
        use_distance_mode: bool = False,
        n_proteins: int = 1,
        defaults: dict | None = None,
    ) -> dict:
        """All placement settings of one protein in the sidebar.

        - **Copies** (only with one protein): how many, and whether they move together on a
          grid or each copy is placed independently.
        - **Position**: x and y from the box centre, or random.
        - **Height**: absolute z, or a gap above the membrane surface.
        - **Rotation**: angles about x, y and z, or random (also per system or per copy).
        - **Independent copies**: their own position, height and rotation, each in its own
          expandable section; they start on the grid of "Move together".

        Args:
            protein_name: The protein's file name.
            boxx: Box length x, in nm (limits the position).
            boxy: Box length y, in nm.
            boxz: Box height z, in nm.
            n_systems: Number of systems (offers "re-randomize for each system").
            protein_only: No membrane, so no "Height above Membrane".
            use_distance_mode: Two proteins placed by their gap; x/y are not asked.
            n_proteins: Number of proteins (copies are offered for a single protein only).
            defaults: Saved values to start from.

        Returns:
            The settings of this protein (`more_copies`, `copy_number`, `copy_mode`,
            `copies`, `randomize_pos`, `cx`, `cy`, `z_method`, `cz`, `distance_to_mem`,
            `randomize_rot`, `rx`, `ry`, `rz`, …).
        """
        d = defaults or {}
        ns = _ns(protein_name)

        st.session_state.setdefault(f"{ns}_more_copies", bool(_default(d, "more_copies", False)))
        st.session_state.setdefault(f"{ns}_copy_number", int(_default(d, "copy_number", 2)))
        st.session_state.setdefault(f"{ns}_randomize_pos", bool(_default(d, "randomize_pos", False)))
        st.session_state.setdefault(f"{ns}_randomize_pos_every", bool(_default(d, "randomize_pos_every", False)))
        st.session_state.setdefault(f"{ns}_z_method", _default(d, "z_method", "Absolute z position"))
        st.session_state.setdefault(f"{ns}_randomize_rot", bool(_default(d, "randomize_rot", False)))
        st.session_state.setdefault(f"{ns}_randomize_rot_every", bool(_default(d, "randomize_rot_every", False)))
        st.session_state.setdefault(f"{ns}_randomize_rot_all_copies", bool(_default(d, "randomize_rot_all_copies", False)))
        st.session_state.setdefault(f"{ns}_copy_mode", _default(d, "copy_mode", "together"))
        # cx, cy, cz, rx, ry, rz, distance_to_mem are seeded inside _num_input

        more_copies = False
        copy_number = 2
        copy_mode = "together"
        if n_proteins == 1:
            st.sidebar.markdown("---")
            more_copies = st.sidebar.checkbox(
                "Have more than 1 copy of the protein inside the system",
                key=f"{ns}_more_copies",
            )
            if more_copies:
                st.sidebar.subheader("Multiple protein setup")
                st.session_state[f"{ns}_copy_number"] = max(2, int(st.session_state[f"{ns}_copy_number"]))
                copy_number = st.sidebar.number_input(
                    "Number of copies per system",
                    min_value=2,
                    step=1,
                    key=f"{ns}_copy_number",
                )
                copy_mode = st.sidebar.radio(
                    "Copies",
                    ["together", "independent"],
                    format_func=lambda m: {"together": "Move together (spread on a grid)",
                                           "independent": "Move each copy independently"}[m],
                    key=f"{ns}_copy_mode",
                )

        st.sidebar.subheader("Protein Placement")
        cx, cy = 0.0, 0.0
        randomize_pos = False
        randomize_pos_every = False

        if not more_copies and not use_distance_mode:
            randomize_pos = st.sidebar.checkbox("Randomize xy position", key=f"{ns}_randomize_pos")
            if randomize_pos and n_systems > 1:
                randomize_pos_every = st.sidebar.checkbox(
                    "Re-randomize position for each system", key=f"{ns}_randomize_pos_every"
                )
            if not randomize_pos:
                cx = _num_input("Position X (nm)", f"{ns}_cx", -boxx / 2, boxx / 2, 0.1, _default(d, "cx", 0.0))
                cy = _num_input("Position Y (nm)", f"{ns}_cy", -boxy / 2, boxy / 2, 0.1, _default(d, "cy", 0.0))
        elif use_distance_mode:
            st.sidebar.info("X/Y position set automatically via protein distance.")

        cz, distance_to_mem = 0.0, 2.0
        if not protein_only:
            z_method = st.sidebar.selectbox(
                "Z placement method",
                ["Absolute z position", "Height above Membrane"],
                key=f"{ns}_z_method",
            )
        else:
            z_method = "Absolute z position"

        independent = more_copies and copy_mode == "independent"
        if z_method == "Absolute z position":
            if not independent:  # independent copies have their own z below
                cz = _num_input("Absolute Z (nm)", f"{ns}_cz", -boxz / 2, boxz / 2, 0.1, _default(d, "cz", 0.0))
        else:
            distance_to_mem = _num_input(
                "Distance above membrane (nm)", f"{ns}_distance_to_mem",
                0.0, boxz / 2, 0.1, _default(d, "distance_to_mem", 2.0),
            )

        # Rotation
        st.sidebar.subheader("Protein Rotation")
        rx, ry, rz = 0.0, 0.0, 0.0
        randomize_rot = st.sidebar.checkbox("Randomize rotation", key=f"{ns}_randomize_rot")
        randomize_rot_every = False
        randomize_rot_all_copies = False
        if randomize_rot and n_systems > 1:
            randomize_rot_every = st.sidebar.checkbox(
                "Re-randomize rotation for each system", key=f"{ns}_randomize_rot_every"
            )
        if more_copies and not independent:
            randomize_rot_all_copies = st.sidebar.checkbox(
                "Randomize the rotation for all copies", key=f"{ns}_randomize_rot_all_copies"
            )
        if not randomize_rot and not independent:
            rx = _num_input("Rotation X (°)", f"{ns}_rx", -180.0, 180.0, 1.0, _default(d, "rx", 0.0))
            ry = _num_input("Rotation Y (°)", f"{ns}_ry", -180.0, 180.0, 1.0, _default(d, "ry", 0.0))
            rz = _num_input("Rotation Z (°)", f"{ns}_rz", -180.0, 180.0, 1.0, _default(d, "rz", 0.0))

        # Independent copies: every copy has its own position (and rotation, unless randomized).
        # They start on the grid used by "Move together", so switching modes keeps the layout.
        copies = []
        if independent:
            from openmembranebuilder.builders.membrane.schematic import make_grid
            saved = d.get("copies") or []
            grid = make_grid(boxx, boxy, int(copy_number))
            st.sidebar.subheader("Copies")
            for k in range(int(copy_number)):
                s0 = saved[k] if k < len(saved) else {}
                gx, gy = grid[k] if k < len(grid) else (0.0, 0.0)
                with st.sidebar.expander(f"Copy {k + 1}", expanded=k == 0):
                    c = {
                        "cx": _num_input("Position X (nm)", f"{ns}_c{k}_cx", -boxx / 2, boxx / 2, 0.1, _default(s0, "cx", gx), container=st),
                        "cy": _num_input("Position Y (nm)", f"{ns}_c{k}_cy", -boxy / 2, boxy / 2, 0.1, _default(s0, "cy", gy), container=st),
                        "cz": 0.0, "rx": 0.0, "ry": 0.0, "rz": 0.0,
                    }
                    if z_method == "Absolute z position":
                        c["cz"] = _num_input("Absolute Z (nm)", f"{ns}_c{k}_cz", -boxz / 2, boxz / 2, 0.1, _default(s0, "cz", 0.0), container=st)
                    if not randomize_rot:
                        for ax in ("rx", "ry", "rz"):
                            c[ax] = _num_input(f"Rotation {ax[1].upper()} (°)", f"{ns}_c{k}_{ax}", -180.0, 180.0, 1.0, _default(s0, ax, 0.0), container=st)
                copies.append(c)

        return {
            "more_copies": more_copies,
            "copy_number": copy_number,
            "copy_mode": copy_mode if more_copies else "together",
            "copies": copies,
            "randomize_pos": randomize_pos,
            "randomize_pos_every": randomize_pos_every,
            "cx": cx,
            "cy": cy,
            "z_method": z_method,
            "cz": cz,
            "distance_to_mem": distance_to_mem,
            "randomize_rot": randomize_rot,
            "randomize_rot_every": randomize_rot_every,
            "randomize_rot_all_copies": randomize_rot_all_copies,
            "rx": rx,
            "ry": ry,
            "rz": rz,
        }
    @staticmethod
    def render_multi_protein_params(
        protein_files: dict,           
        boxx, boxy, boxz, n_systems, protein_only
    ) -> dict:
        """Placement settings for all proteins (earlier version, replaced by the one below)."""
        if not protein_files:
            return {}  
        results = {}
        use_distance = False

        if len(protein_files) == 2:
            st.sidebar.markdown("---")
            use_distance = st.sidebar.checkbox(
                "Set distance between proteins instead of X/Y positions",
                key="use_protein_distance"
            )
            if use_distance:
                st.sidebar.number_input(
                    "Distance between proteins (nm)", 0.0, min(boxx, boxy),
                    value=st.session_state.get("protein_distance", 1.0),
                    step=0.1, key="protein_distance"
                )

        for pdb_name in protein_files:
            st.sidebar.markdown("---")
            st.sidebar.markdown("---")
            st.sidebar.markdown(f"## {Path(pdb_name).stem} -- Params")
            results[pdb_name] = ProteinGui.render_insertion_params_for_protein(
                protein_name=pdb_name,
                boxx=boxx, boxy=boxy, boxz=boxz,
                n_systems=n_systems,
                protein_only=protein_only,
                use_distance_mode=(use_distance and len(protein_files) == 2),
                n_proteins=len(protein_files),   
            )
        st.session_state["insertion_params"] = results
        return results    

    @staticmethod
    def render_multi_protein_params(
        protein_files: dict, boxx, boxy, boxz, n_systems, protein_only, config=None
    ) -> dict:
        """Placement settings for all proteins, one block per protein.

        With exactly two proteins you can also place them by the gap between them instead of
        by x/y positions.

        Args:
            protein_files: The proteins.
            boxx: Box length x, in nm.
            boxy: Box length y, in nm.
            boxz: Box height z, in nm.
            n_systems: Number of systems.
            protein_only: Proteins in water, without a membrane.
            config: The saved settings, used to fill the widgets.

        Returns:
            The settings of every protein, by file name.
        """
        if not protein_files:
            return {}

        # saved values from the config are the defaults of the widgets
        saved_poses = (getattr(config, "protein_params", None) or {}).get("R0001", {})
        saved_ip = getattr(config, "insertion_params", None) or {}

        results = {}

        use_distance = False
        if len(protein_files) == 2:
            st.session_state.setdefault(
                "use_protein_distance", bool(getattr(config, "use_protein_distance", False))
            )
            st.session_state.setdefault(
                "protein_distance", float(getattr(config, "protein_distance", 1.0))
            )
            max_dist = float(min(boxx, boxy))
            st.session_state["protein_distance"] = min(max(float(st.session_state["protein_distance"]), 0.0), max_dist)

            st.sidebar.markdown("---")
            use_distance = st.sidebar.checkbox(
                "Set distance between proteins instead of X/Y positions",
                key="use_protein_distance",
            )
            if use_distance:
                st.sidebar.number_input(
                    "Distance between proteins (nm)",
                    min_value=0.0,
                    max_value=max_dist,
                    step=0.1,
                    key="protein_distance",
                )

        for pdb_name in protein_files:
            st.sidebar.markdown("---")
            st.sidebar.markdown(f"## {Path(pdb_name).stem} -- Params")
            defaults = {**saved_poses.get(pdb_name, {}), **saved_ip.get(pdb_name, {})}
            results[pdb_name] = ProteinGui.render_insertion_params_for_protein(
                protein_name=pdb_name,
                boxx=boxx,
                boxy=boxy,
                boxz=boxz,
                n_systems=n_systems,
                protein_only=protein_only,
                use_distance_mode=(use_distance and len(protein_files) == 2),
                n_proteins=len(protein_files),
                defaults=defaults,
            )

        st.session_state["insertion_params"] = results
        return results


# ─── Page renderers ───────────────────────────────────────────────────────────

def render_membrane(force_fields, ff_data, global_info, config, temp_folder, protein_only: bool = False) -> dict | str | None:
    """Draw the membrane page.

    The first time the page is entered, its widgets are filled from the saved settings.

    Args:
        force_fields: Names of the available force fields.
        ff_data: Lipids and lipid descriptions of every force field.
        global_info: The standard folder locations.
        config: The saved settings.
        temp_folder: Folder for uploaded templates.
        protein_only: Proteins in water: the composition part is not shown.
    """
    if not st.session_state.get("_membrane_seeded"):
        st.session_state["selected_ff"] = config.selected_ff
        st.session_state["box_type"] = config.box_type
        st.session_state["boxx"] = config.box_x
        st.session_state["boxy"] = config.box_y
        st.session_state["boxz"] = config.box_z
        st.session_state["n_systems"] = config.n_systems
        st.session_state["output_name"] = config.output_name
        st.session_state["abs_lip_vals"] = config.abs_lip_vals
        # entries may be a dict {template: entries}; the manual editor needs a list
        if isinstance(config.entries, list) and config.entries:
            st.session_state["entries"] = config.entries
            _reset_entry_widgets()
        solv = _parse_solvation(config.solvation)
        st.session_state["pos_ion"] = solv["pos_ion"]
        st.session_state["neg_ion"] = solv["neg_ion"]
        st.session_state["salt_molarity"] = solv["salt_molarity"]
        st.session_state["replace_w_percentage"] = solv["replace_w_percentage"]
        st.session_state["_membrane_seeded"] = True

    result = {}

    def sidebar():
        """Sidebar: force field, box, solvation, number of systems, output name."""
        result["selected_ff"] = MembraneGui.render_force_field_selector(force_fields)
        if not protein_only or st.session_state.get("protein_files"):
            _ff_msg = ff_mismatch_message(st.session_state.get("protein_files", {}), result["selected_ff"])
            if _ff_msg:
                st.sidebar.warning(_ff_msg)
        (
            result["selected_box_type"],
            result["box_x"],
            result["box_y"],
            result["box_z"],
        ) = MembraneGui.render_box_params()
        result["solvation"] = MembraneGui.render_solvation_input(result["selected_ff"])
        result["n_systems"] = MembraneGui.render_n_systems_input()
        result["output_name"] = MembraneGui.render_output_name_input()

    def content():
        """Main area: errors, template choice, lipid table, Back / Next."""
        for err in st.session_state.get("_membrane_errors", []):
            st.error(err)

        selected_ff = result.get("selected_ff", force_fields[0])
        available_lipids = ff_data[selected_ff]["available_lipids"]
        lipid_map = ff_data[selected_ff]["lipid_map"]

        if not protein_only:
            template_source = st.session_state.get("template_source", {
                "kind": "manual",
                "name": "",
                "path": "",
                "loaded_for": "",
                "entries": [],
                "templates": {},
            })
            result["template_active"] = MembraneGui.render_template_uploader(
                available_lipids,
                global_info,
                temp_folder,
                template_source=template_source,
                entries=st.session_state.get("entries", []),
                abs_lip_vals=st.session_state.get("abs_lip_vals", False),
            )
            from openmembranebuilder.builders.membrane.page import _hydrate_template_source
            template_source = _hydrate_template_source(available_lipids, global_info)
            st.session_state["template_source"] = template_source
            result["template_path"] = template_source.get("path", "")
            MembraneGui.render_lipid_mapping(lipid_map)
            result["entries"] = MembraneGui.render_membrane_entries(
                available_lipids,
                entries=st.session_state.get("entries", []),
                abs_lip_vals=st.session_state.get("abs_lip_vals", False),
                template_source=template_source,
            )
            result["abs_lip_vals"] = st.session_state.get("abs_lip_vals", False)

        with st.container(key="nav_membrane"):
            st.divider()
            left, _, right = st.columns([1, 6, 1])
            if left.button("← Back", use_container_width=True):
                st.session_state["_membrane_action"] = "back"
                st.rerun()
            if right.button("Next →", use_container_width=True, type="primary"):
                st.session_state["_membrane_action"] = "next"
                st.rerun()

    GlobalGui.render_page(title="Membrane Builder", content_fn=content, sidebar_fn=sidebar, layout="wide")

    # (never reached: render_page ends the page drawing)
    GlobalGui.render_page(title="Membrane Builder", content_fn=content,sidebar_fn=sidebar, layout="wide")

def render_protein_insertion(global_info, boxx, boxy, boxz, n_systems, config, temp_folder, protein_only: bool = False) -> dict | str | None:
    """Draw the protein insertion page.

    The first time the page is entered, its widgets are filled from the saved settings
    (poses, options, independent copies, box).

    Args:
        global_info: The standard folder locations.
        boxx: Box length x, in nm.
        boxy: Box length y, in nm.
        boxz: Box height z, in nm.
        n_systems: Number of systems.
        config: The saved settings.
        temp_folder: Folder for uploads.
        protein_only: Proteins in water, without a membrane.
    """
    if not st.session_state.get("_protein_seeded"):
        for pdb_name, pose in config.protein_params.get("R0001", {}).items():
            ns = _ns(pdb_name)
            st.session_state[f"{ns}_cx"] = pose.get("cx", 0.0)
            st.session_state[f"{ns}_cy"] = pose.get("cy", 0.0)
            st.session_state[f"{ns}_cz"] = pose.get("cz", 0.0)
            st.session_state[f"{ns}_rx"] = pose.get("rx", 0.0)
            st.session_state[f"{ns}_ry"] = pose.get("ry", 0.0)
            st.session_state[f"{ns}_rz"] = pose.get("rz", 0.0)

        saved_ip = getattr(config, "insertion_params", {}) or {}
        for pdb_name, ip in saved_ip.items():
            ns = _ns(pdb_name)
            st.session_state[f"{ns}_randomize_rot"] = ip.get("randomize_rot", False)
            st.session_state[f"{ns}_randomize_rot_every"] = ip.get("randomize_rot_every", False)
            st.session_state[f"{ns}_randomize_rot_all_copies"] = ip.get("randomize_rot_all_copies", False)
            st.session_state[f"{ns}_randomize_pos"] = ip.get("randomize_pos", False)
            st.session_state[f"{ns}_randomize_pos_every"] = ip.get("randomize_pos_every", False)
            st.session_state[f"{ns}_z_method"] = ip.get("z_method", "Absolute z position")
            st.session_state[f"{ns}_distance_to_mem"] = _default(ip, "distance_to_mem", 2.0)
            st.session_state[f"{ns}_more_copies"] = ip.get("more_copies", False)
            st.session_state[f"{ns}_copy_number"] = ip.get("copy_number", 2)
            st.session_state[f"{ns}_copy_mode"] = ip.get("copy_mode", "together")
            for k, cp in enumerate(ip.get("copies") or []):
                for ax, v in cp.items():
                    st.session_state[f"{ns}_c{k}_{ax}"] = v

        st.session_state["use_protein_distance"] = getattr(config, "use_protein_distance", False)
        st.session_state["protein_distance"] = getattr(config, "protein_distance", 1.0)
        # the box can be adjusted here too (same widgets as on the membrane page)
        st.session_state["box_type"] = config.box_type
        st.session_state["boxx"] = config.box_x
        st.session_state["boxy"] = config.box_y
        st.session_state["boxz"] = config.box_z
        st.session_state["_protein_seeded"] = True

    result = {}

    def sidebar():
        """Sidebar: box, list of proteins, placement settings of each protein."""
        result["box_type"], result["boxx"], result["boxy"], result["boxz"] = MembraneGui.render_box_params()
        result["protein_files"] = ProteinGui.render_protein_list()
        protein_files = result["protein_files"]
        result["insertion_params"] = ProteinGui.render_multi_protein_params(
            protein_files,
            result["boxx"],
            result["boxy"],
            result["boxz"],
            n_systems,
            protein_only=protein_only,
            config=config,
        )

    def content():
        """Main area: system view, Back / Next."""
        render_system_view(result, protein_only)
        with st.container(key="nav_protein"):
            st.divider()
            left, _, right = st.columns([1, 6, 1])
            if left.button("← Back", use_container_width=True):
                st.session_state["_protein_action"] = "back"
                st.rerun()
            if right.button("Next →", use_container_width=True, type="primary"):
                protein_files = st.session_state.get("protein_files", {})
                if not protein_files:
                    st.error("No proteins available. Go back to the martinize builder and add at least one protein.")
                else:
                    st.session_state["_protein_action"] = "next"
                    st.rerun()

    GlobalGui.render_page(
        title="Protein Insertion",
        content_fn=content,
        sidebar_fn=sidebar,
        layout="wide",
    )


@st.cache_data(show_spinner=False)
def _bead_coords(pdb_path: str, mtime: float):
    """The bead positions of a protein file, in Å.

    Remembered per file and modification time, so the file is read only once.

    Args:
        pdb_path: The protein structure.
        mtime: The file's modification time (a changed file is read again).

    Returns:
        The positions, one row per bead.
    """
    import MDAnalysis as mda
    return mda.Universe(pdb_path).atoms.positions.astype(float)


def render_system_view(result: dict, protein_only: bool) -> None:
    """The schematic system view of the protein insertion page, with its warnings.

    Shows box, membrane and every protein (or copy) as COBY will place them, seen from
    the side (x–z, y–z) or from the top (x–y). Warns when a protein sticks out of the box,
    and when proteins or copies overlap (beads closer than 0.3 nm).

    Args:
        result: The values of the sidebar (box and placement settings).
        protein_only: Proteins in water: no membrane is drawn.
    """
    import numpy as np
    from openmembranebuilder.builders.membrane import schematic as sch

    box_type = result.get("box_type", "rectangular")
    extent = sch.box_extent(box_type, result["boxx"], result["boxy"], result["boxz"])
    protein_files = result.get("protein_files", {}) or {}
    params = result.get("insertion_params", {}) or {}
    use_distance = bool(st.session_state.get("use_protein_distance")) and len(protein_files) == 2
    distance = float(st.session_state.get("protein_distance", 1.0))

    proteins, coords = [], {}
    for name, paths in protein_files.items():
        try:
            coords[name] = _bead_coords(paths["pdb_path"], os.path.getmtime(paths["pdb_path"]))
        except Exception as exc:  # noqa: BLE001 - a protein that cannot be read is simply not drawn
            st.caption(f"{name}: cannot draw ({exc})")
    for k, (name, c) in enumerate(coords.items()):
        p = params.get(name, {})
        rx, ry, rz = (p.get("rx", 0.0), p.get("ry", 0.0), p.get("rz", 0.0))
        cx, cy = p.get("cx", 0.0), p.get("cy", 0.0)
        if use_distance:  # as the runner: surfaces `distance` apart along x, first protein at +x
            cx = (1 if k == 0 else -1) * (sch.half_extent_x(c, rx, ry, rz) + distance / 2)
            cy = 0.0
        if p.get("z_method") == "Height above Membrane" and not protein_only:
            cz = sch.height_above_membrane_cz(c, rx, ry, rz, p.get("distance_to_mem", 2.0))
        else:
            cz = p.get("cz", 0.0)
        independent = p.get("more_copies") and p.get("copy_mode") == "independent"
        if independent:
            placed = []
            for cp in p.get("copies", []):
                crx, cry, crz = cp.get("rx", 0.0), cp.get("ry", 0.0), cp.get("rz", 0.0)
                ccz = (sch.height_above_membrane_cz(c, crx, cry, crz, p.get("distance_to_mem", 2.0))
                       if p.get("z_method") == "Height above Membrane" and not protein_only else cp.get("cz", 0.0))
                placed.append(sch.place(c, cp.get("cx", 0.0), cp.get("cy", 0.0), ccz, crx, cry, crz))
        elif p.get("more_copies"):
            n = int(p.get("copy_number", 2))
            placed = [sch.place(c, x, y, cz, rx, ry, rz) for x, y in sch.make_grid(extent[0], extent[1], n)[:n]]
        else:
            placed = [sch.place(c, cx, cy, cz, rx, ry, rz)]
        proteins.append({
            "name": os.path.splitext(name)[0],
            "copies": placed,
            # dashed = only one possible arrangement (random, or grid cells picked at random)
            "dashed": bool(p.get("randomize_pos") or p.get("randomize_rot") or (p.get("more_copies") and not independent)),
        })

    st.subheader("System view")
    view = st.radio("View", list(sch.VIEWS), horizontal=True, key="protein_view", label_visibility="collapsed")
    st.markdown(sch.system_view_svg(extent, proteins, membrane=not protein_only, view=view), unsafe_allow_html=True)
    notes = ["Proteins are drawn as the outline of their beads, centred and rotated as COBY places them."]
    if any(p["dashed"] for p in proteins):
        notes.append("Dashed: random position/rotation or copies (one possible arrangement shown).")
    if not protein_only:
        notes.append("Membrane drawn with a typical Martini thickness; 'Height above Membrane' uses the real membrane at build time.")
    if box_type != "rectangular":
        notes.append("Hexagonal boxes are drawn as their bounding box.")
    st.caption(" ".join(notes))
    half = np.array(extent) / 2
    out = [p["name"] for p in proteins if any((np.abs(c) > half).any() for c in p["copies"])]
    if out:
        st.warning(f"{', '.join(out)} stick(s) out of the box: enlarge the box or move/rotate the protein.")
    # overlapping proteins or copies (beads closer than 0.3 nm)
    from scipy.spatial import cKDTree
    bodies = [(f"{p['name']} {k + 1}" if len(p["copies"]) > 1 else p["name"], c) for p in proteins for k, c in enumerate(p["copies"])]
    clashes = [f"{a} / {b}" for i, (a, ca) in enumerate(bodies) for (b, cb) in bodies[i + 1:]
               if cKDTree(ca).query(cb, distance_upper_bound=0.3)[0].min() < 0.3]
    if clashes and not any(p["dashed"] for p in proteins):
        st.warning("Overlapping: " + ", ".join(clashes) + ". Move or rotate them apart.")


def render_preview(config, has_protein: bool, protein_only: bool = False) -> None:
    """Draw the preview page: all settings at a glance and the Build button.

    When proteins and membrane belong to different Martini generations, a warning is
    shown and building needs an explicit confirmation.

    Args:
        config: The settings to show.
        has_protein: Show the protein part.
        protein_only: Proteins in water: no composition is shown.
    """
    def md(label, val):
        # values are user input (output name, file names) -> escape before injecting HTML
        """One label with its value below it (text is made safe for the page)."""
        return (
            f'<p style="font-size:17px; color:gray; margin:0">{html.escape(str(label))}</p>'
            f'<p style="font-size:20px; font-weight:600; margin:0 0 16px 0">'
            f'{html.escape(str(val)) if val else "&nbsp;"}</p>'
        )

    def row(label1, val1, label2="", val2=""):
        """Two label/value pairs side by side."""
        c1, c2 = st.columns(2)
        c1.markdown(md(label1, val1), unsafe_allow_html=True)
        c2.markdown(md(label2, val2), unsafe_allow_html=True)

    def row4(l1, v1, l2, v2, l3, v3, l4, v4):
        """Four label/value pairs side by side."""
        c1, c2, c3, c4 = st.columns(4)
        for c, l, v in [(c1, l1, v1), (c2, l2, v2), (c3, l3, v3), (c4, l4, v4)]:
            c.markdown(md(l, v), unsafe_allow_html=True)

    def content():
        """The summary sections and the Back / Build buttons."""
        st.subheader("Please check your settings before building.")

        with st.expander("⚙️ General", expanded=True):
            row("Force field", config.selected_ff, "Number of systems", str(config.n_systems))
            if config.output_name:
                row("Output name", config.output_name)

        with st.expander("📦 Box", expanded=True):
            row4("Type", config.box_type, "X (nm)", str(config.box_x), "Y (nm)", str(config.box_y), "Z (nm)", str(config.box_z))

        with st.expander("💧 Solvation", expanded=True):
            st.code(config.solvation or "—")

        if not protein_only:
            with st.expander("🧪 Membrane Composition", expanded=True):
                row("Absolute values", str(config.abs_lip_vals), "Template active", str(config.template_active))
                headers = ["Lipid", "NU" if config.abs_lip_vals else "RU", "NL" if config.abs_lip_vals else "RL", "AU", "AL"]
                if config.template_active and isinstance(config.entries, dict):
                    for template_name, entries in config.entries.items():
                        st.markdown(f"**{template_name}**")
                        st.dataframe(pd.DataFrame(entries, columns=headers), use_container_width=True, hide_index=True)
                elif config.entries:
                    st.dataframe(pd.DataFrame(config.entries, columns=headers), use_container_width=True, hide_index=True)

        if has_protein:
            with st.expander("🧬 Protein Insertion", expanded=True):
                distance_mode = getattr(config, "use_protein_distance", False) and len(config.protein_files) == 2
                if distance_mode:
                    row("Distance between proteins (nm)", str(config.protein_distance))
                    st.divider()

                insertion_params = getattr(config, "insertion_params", {}) or {}
                for pdb_name, paths in config.protein_files.items():
                    p = config.protein_params.get("R0001", {}).get(pdb_name, {})
                    ip = insertion_params.get(pdb_name, {})

                    st.markdown(f"#### {pdb_name}")
                    row("PDB file", pdb_name, "ITP file", os.path.basename(paths["itp_path"]))

                    if ip.get("more_copies") and ip.get("copy_number", 1) > 1:
                        row(
                            "Copies per system",
                            str(ip["copy_number"]),
                            "Randomize rotation (all copies)",
                            str(ip.get("randomize_rot_all_copies", False)),
                        )

                    if distance_mode:
                        pass
                    elif ip.get("randomize_pos"):
                        row(
                            "Position",
                            "Randomized",
                            "Re-randomize per system" if config.n_systems > 1 else "",
                            str(ip.get("randomize_pos_every", False)) if config.n_systems > 1 else "",
                        )
                    else:
                        row("Position X, Y (nm)", f"{p.get('cx', 0.0)}, {p.get('cy', 0.0)}")

                    if ip.get("z_method") == "Height above Membrane":
                        row(
                            "Z placement",
                            "Height above Membrane",
                            "Distance to membrane (nm)",
                            str(ip.get("distance_to_mem", 2.0)),
                        )
                    else:
                        row("Z placement", "Absolute", "Z (nm)", str(p.get("cz", 0.0)))

                    if ip.get("randomize_rot"):
                        row(
                            "Rotation",
                            "Randomized",
                            "Re-randomize per system" if config.n_systems > 1 else "",
                            str(ip.get("randomize_rot_every", False)) if config.n_systems > 1 else "",
                        )
                    else:
                        row("Rotation X, Y, Z (°)", f"{p.get('rx', 0.0)}, {p.get('ry', 0.0)}, {p.get('rz', 0.0)}")

                    st.divider()

        # Martini 2 / 3 mismatch: warn, and require an explicit confirmation before building
        blocked = False
        ff_msg = ff_mismatch_message(config.protein_files if has_protein else {}, config.selected_ff)
        if ff_msg:
            if st.session_state.get("_ff_mismatch_msg") != ff_msg:
                st.session_state["_ff_mismatch_msg"] = ff_msg
                st.session_state["ff_mismatch_ack"] = False
            st.warning(ff_msg)
            ack = st.checkbox("I understand the risk and want to build anyway", key="ff_mismatch_ack")
            blocked = not ack
        else:
            st.session_state.pop("_ff_mismatch_msg", None)

        with st.container(key="nav_preview"):
            st.divider()
            left, _, right = st.columns([1, 6, 1])
            if left.button("← Back", use_container_width=True):
                st.session_state["_preview_action"] = "back"
                st.rerun()
            if right.button("🔨 Build", use_container_width=True, type="primary", disabled=blocked):
                st.session_state["_preview_action"] = "build"
                st.rerun()

    GlobalGui.render_page(title="Preview", content_fn=content, layout="wide")


def _parse_solvation(solvation: str) -> dict:
    """Read ions, salt and antifreeze water back from a solvation text.

    Args:
        solvation: Solvation in COBY's notation.

    Returns:
        `pos_ion`, `neg_ion`, `salt_molarity` and `replace_w_percentage`, with defaults
        for anything not found.
    """
    result = {
        "pos_ion": "NA",
        "neg_ion": "CL",
        "salt_molarity": 0.15,
        "replace_w_percentage": 0.0,
    }
    if not solvation:
        return result

    pos = re.search(r"pos:(\w+)", solvation)
    neg = re.search(r"neg:(\w+)", solvation)
    salt = re.search(r"salt_molarity:([\d.]+)", solvation)
    wf = re.search(r"solv:WF([\d.]+)", solvation)

    if pos:
        result["pos_ion"] = pos.group(1)
    if neg:
        result["neg_ion"] = neg.group(1)
    if salt:
        result["salt_molarity"] = float(salt.group(1))
    if wf:
        result["replace_w_percentage"] = float(wf.group(1))
    return result
