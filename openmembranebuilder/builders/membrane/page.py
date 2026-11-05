import copy
import os
from datetime import datetime
from enum import Enum, auto

import streamlit as st

from openmembranebuilder.builders.membrane.gui import _reset_entry_widgets, render_membrane, render_protein_insertion, render_preview
from openmembranebuilder.builders.membrane.parser import MembraneParser
from openmembranebuilder.builders.membrane.config import Config
from openmembranebuilder.core.global_info import GlobalInfo
from openmembranebuilder.builders.membrane.runner import MembraneRunner


class MembraneStep(Enum):
    MEMBRANE = auto()
    PROTEIN_INSERTION = auto()
    PREVIEW = auto()
    RESULTS = auto()


def _get_step() -> MembraneStep:
    if "membrane_step" not in st.session_state:
        st.session_state.membrane_step = MembraneStep.MEMBRANE
    return st.session_state.membrane_step


def _set_step(step: MembraneStep):
    st.session_state.membrane_step = step
    st.session_state.pop("_membrane_errors", None)
    if step == MembraneStep.MEMBRANE:
        st.session_state.pop("_membrane_seeded", None)
    elif step == MembraneStep.PROTEIN_INSERTION:
        st.session_state.pop("_protein_seeded", None)


def main(
    global_info,
    has_protein=False,
    temp_folder="",
    protein_only=False,
    extern_config=None,
):
    extern_config = extern_config or {}

    parser = MembraneParser()

    extern_signature = (
        tuple(
            sorted(
                (name, tuple(sorted(info.items())))
                for name, info in extern_config.get("lipid_definitions", {}).items()
            )
        ),
        tuple(sorted(extern_config.get("molecule_builder_index", {}).items())),
        tuple(sorted(extern_config.get("config_overrides", {}).items())),
    )

    if (
        "ff_data" not in st.session_state
        or "force_fields" not in st.session_state
        or st.session_state.get("_ff_data_extern_signature") != extern_signature
    ):
        force_fields = parser.get_force_fields()
        st.session_state["ff_data"] = {
            ff: {
                "available_lipids": parser.get_available_lipids(
                    ff,
                    extern_config=extern_config,
                ),
                "lipid_map": parser.build_lipid_map(
                    ff,
                    extern_config=extern_config,
                ),
            }
            for ff in force_fields
        }
        st.session_state["force_fields"] = force_fields
        st.session_state["_ff_data_extern_signature"] = extern_signature

    if "config" not in st.session_state:
        st.session_state["config"] = Config(base_folder=global_info.base_folder)
    st.session_state["has_protein"]  = has_protein
    st.session_state["protein_only"] = protein_only

    st.session_state["has_protein"] = has_protein
    st.session_state["protein_only"] = protein_only

    ff_data = st.session_state["ff_data"]
    force_fields = st.session_state["force_fields"]
    config = st.session_state["config"]

    _apply_extern_config(config, extern_config)

    for rep_key, rep_val in config.protein_params.items():
        if rep_val and not isinstance(next(iter(rep_val.values())), dict):
            config.protein_params[rep_key] = {}

    step = _get_step()

    if step == MembraneStep.MEMBRANE:
        result = _render_membrane(
            global_info,
            ff_data,
            force_fields,
            has_protein,
            config,
            temp_folder,
            protein_only=protein_only,
        )
        if result == "back":
            return "back", None

    elif step == MembraneStep.PROTEIN_INSERTION:
        _render_protein_insertion(
            global_info,
            config,
            temp_folder,
            protein_only=protein_only,
        )

    elif step == MembraneStep.PREVIEW:
        action = st.session_state.pop("_preview_action", None)

        if action == "back":
            back_step = MembraneStep.PROTEIN_INSERTION if has_protein else MembraneStep.MEMBRANE
            _set_step(back_step)
            st.rerun()

        if action == "build":
            st.session_state["_building"] = True
            _set_step(MembraneStep.RESULTS)
            st.rerun()

        render_preview(config, has_protein, protein_only=protein_only)

    elif step == MembraneStep.RESULTS:
        if st.session_state.get("_building"):
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            folder_name = (
                f"OpenMembraneBuilder-{timestamp}-{config.output_name}"
                if config.output_name else f"OpenMembraneBuilder-{timestamp}"
            )
            config.base_folder = os.path.join(global_info.base_folder, folder_name)
            config.output_name  = config.base_folder

            with st.spinner("Building membrane system..."):
                output_pdbs, zip_path = MembraneRunner().run(config)

            st.session_state["output_pdbs"] = output_pdbs
            st.session_state["zip_path"] = zip_path
            st.session_state["config"] = config
            st.session_state.pop("_building", None)

        return "done", config

    return None, None


def _apply_extern_config(config, extern_config: dict) -> None:
    overrides = extern_config.get("config_overrides", {})
    for key, value in overrides.items():
        setattr(config, key, value)

    config.extern_config = extern_config
    config.molecule_builder_index = extern_config.get("molecule_builder_index", {})


def _default_template_source() -> dict:
    return {
        "kind": "manual",
        "name": "",
        "path": "",
        "loaded_for": "",
        "entries": [],
        "templates": {},
        "manual_backup": [],
        "abs_lip_vals": False,
    }


def _get_template_source() -> dict:
    state = st.session_state.get("template_source")
    if isinstance(state, dict):
        return state
    default = _default_template_source()
    st.session_state["template_source"] = default
    return default


def _set_template_source(template_source: dict) -> dict:
    st.session_state["template_source"] = template_source
    return template_source


def _save_manual_backup(entries=None) -> dict:
    template_source = _get_template_source()
    current_entries = entries if isinstance(entries, list) else st.session_state.get("entries")
    template_source["manual_backup"] = copy.deepcopy(current_entries) if isinstance(current_entries, list) else []
    return _set_template_source(template_source)


def _restore_manual_backup() -> dict:
    template_source = _get_template_source()
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
    st.session_state["abs_lip_vals"] = bool(template_source.get("abs_lip_vals", False))
    return _set_template_source(template_source)


def _apply_template_selection(kind: str, name: str, path: str, entries: list | None = None, abs_lip_vals: bool | None = None) -> dict:
    template_source = _get_template_source()
    template_source.update({
        "kind": kind,
        "name": name,
        "path": path,
        "entries": list(entries) if isinstance(entries, list) else [],
        "templates": {},
        "loaded_for": path if kind in {"builtin", "upload"} and path else "",
        "abs_lip_vals": bool(abs_lip_vals) if abs_lip_vals is not None else bool(template_source.get("abs_lip_vals", False)),
    })
    if isinstance(entries, list):
        st.session_state["entries"] = list(entries)
    if abs_lip_vals is not None:
        st.session_state["abs_lip_vals"] = bool(abs_lip_vals)
    return _set_template_source(template_source)


def _hydrate_template_source(available_lipids: list, global_info=None) -> dict:
    template_source = _get_template_source()
    kind = template_source.get("kind", "manual")

    if kind not in {"builtin", "upload"}:
        template_source["entries"] = st.session_state.get("entries", [])
        template_source["templates"] = {}
        template_source["loaded_for"] = ""
        template_source["manual_backup"] = template_source.get("manual_backup", [])
        return _set_template_source(template_source)

    path = template_source.get("path", "")
    if not path:
        template_source["entries"] = []
        template_source["loaded_for"] = ""
        return _set_template_source(template_source)

    current_entries = st.session_state.get("entries")
    if template_source.get("loaded_for") == path and template_source.get("entries"):
        return template_source
    if kind in {"builtin", "upload"} and path and isinstance(current_entries, list) and current_entries:
        template_source["entries"] = list(current_entries)
        template_source["loaded_for"] = path
        return _set_template_source(template_source)

    parser = MembraneParser()
    template_entries = parser.load_membrane_template_from_csv(path, available_lipids, gui=False)
    if not template_entries:
        return template_source

    template_source["templates"] = template_entries
    template_source["entries"] = list(next(iter(template_entries.values())))
    template_source["loaded_for"] = path
    _reset_entry_widgets()
    st.session_state["entries"] = list(template_source["entries"])
    st.session_state["abs_lip_vals"] = bool(template_source.get("abs_lip_vals", st.session_state.get("abs_lip_vals", False)))
    return _set_template_source(template_source)


def _render_membrane(global_info, ff_data, force_fields, has_protein, config, temp_folder, protein_only: bool = False) -> str | None:
    action = st.session_state.pop("_membrane_action", None)

    
    # starts. Always copy their session-state values on the action rerun so the
    # values entered by the user, rather than Config defaults, reach the runner.
    if action in ("next", "back"):
        config.output_name = st.session_state.get("output_name", config.output_name)
        config.selected_ff = st.session_state.get("selected_ff", config.selected_ff)
        config.box_x = float(st.session_state.get("boxx", config.box_x))
        config.box_y = float(st.session_state.get("boxy", config.box_y))
        config.box_z = float(st.session_state.get("boxz", config.box_z))
        config.box_type = st.session_state.get("box_type", config.box_type)
        config.n_systems = int(st.session_state.get("n_systems", config.n_systems))
        pos_ion = st.session_state.get("pos_ion", "NA")
        neg_ion = st.session_state.get("neg_ion", "CL")
        salt = st.session_state.get("salt_molarity", 0.15)
        wf = st.session_state.get("replace_w_percentage", 0.0)
        if wf != 0.0:
            config.solvation = (
                f"solv:W{100 - wf} solv:WF{wf} "
                f"params:IMPORTED pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt}"
            )
        else:
            config.solvation = f"solv:W pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt}"

        if protein_only:
            config.entries = []
            config.abs_lip_vals = False
            config.template_active = False
            config.template_path = ""
        else:
            template_source = _hydrate_template_source(ff_data.get(config.selected_ff, {}).get("available_lipids", []), global_info)
            config.abs_lip_vals = st.session_state.get("abs_lip_vals", False)
            config.template_active = template_source.get("kind") in {"builtin", "upload"}
            config.template_path = template_source.get("path", "")
            current_entries = st.session_state.get("entries", template_source.get("entries", []))
            config.entries = current_entries if config.template_active else template_source.get("entries", current_entries)
            st.session_state["entries"] = config.entries
            st.session_state["template_source"] = template_source

    if action == "back":
        st.session_state.pop("_membrane_errors", None)
        return "back"

    if action == "next":
        errors = []
        if not protein_only:
            _, errors = _validate(
                {
                    "entries": config.entries,
                    "abs_lip_vals": config.abs_lip_vals,
                    "template_active": bool(getattr(config, "template_active", False)),
                }
            )

        if errors:
            # Validation failed: stay on this step, remember the errors
            st.session_state["_membrane_errors"] = errors
        else:
            # _set_step also clears any old errors
            _set_step(MembraneStep.PROTEIN_INSERTION if has_protein else MembraneStep.PREVIEW)
            st.rerun()

    
    
    st.session_state.setdefault("template_source", {
        "kind": "manual",
        "name": "",
        "path": "",
        "loaded_for": "",
        "entries": [],
        "templates": {},
    })
    if "entries" not in st.session_state or not isinstance(st.session_state.get("entries"), list):
        st.session_state["entries"] = list(config.entries or [])
    if "abs_lip_vals" not in st.session_state:
        st.session_state["abs_lip_vals"] = bool(config.abs_lip_vals)
    template_source = _get_template_source()
    if isinstance(template_source.get("manual_backup"), list) and not template_source.get("manual_backup") and isinstance(st.session_state.get("entries"), list):
        template_source["manual_backup"] = copy.deepcopy(st.session_state["entries"])
        _set_template_source(template_source)

    render_membrane(force_fields, ff_data, global_info, config, temp_folder, protein_only=protein_only)
    return None

    render_membrane(force_fields, ff_data, global_info, config, temp_folder, protein_only=protein_only)
    return None

def _render_protein_insertion(global_info, config, temp_folder, protein_only: bool = False) -> None:
    action = st.session_state.pop("_protein_action", None)

    if action in ("next", "back"):
        # the box can be changed on this page as well
        config.box_type = st.session_state.get("box_type", config.box_type)
        config.box_x = float(st.session_state.get("boxx", config.box_x))
        config.box_y = float(st.session_state.get("boxy", config.box_y))
        config.box_z = float(st.session_state.get("boxz", config.box_z))
        config.protein_files = st.session_state.get("protein_files", {})
        config.insertion_params = st.session_state.get("insertion_params", {})
        config.use_protein_distance = st.session_state.get("use_protein_distance", False)
        config.protein_distance = st.session_state.get("protein_distance", 1.0)

        # Read per-protein params from insertion_params instead of flat session state
        config.protein_params["R0001"] = {}
        for pdb_name in config.protein_files:
            ns = pdb_name.replace(".", "_").replace(" ", "_")
            ip = config.insertion_params.get(pdb_name, {})
            config.protein_params["R0001"][pdb_name] = {
                "cx": st.session_state.get(f"{ns}_cx", 0.0),
                "cy": st.session_state.get(f"{ns}_cy", 0.0),
                "cz": st.session_state.get(f"{ns}_cz", 0.0),
                "rx": st.session_state.get(f"{ns}_rx", 0.0),
                "ry": st.session_state.get(f"{ns}_ry", 0.0),
                "rz": st.session_state.get(f"{ns}_rz", 0.0),
            }

            # Store per-protein flags onto config for runner and preview
            config.z_method = ip.get("z_method", "Absolute z position")
            config.distance_to_mem = ip.get("distance_to_mem", 2.0)
            config.randomize_pos = ip.get("randomize_pos", False)
            config.randomize_pos_every = ip.get("randomize_pos_every", False)
            config.randomize_rot = ip.get("randomize_rot", False)
            config.randomize_rot_every = ip.get("randomize_rot_every", False)
            config.more_copies = ip.get("more_copies", False)
            config.copy_number = ip.get("copy_number", 1)
            config.copy_mode = ip.get("copy_mode", "together")
            config.copy_params = ip.get("copies", [])
            config.randomize_rot_all_copies = ip.get("randomize_rot_all_copies", False)

    if action == "back":
        _set_step(MembraneStep.MEMBRANE)
        st.rerun()
    if action == "next":
        _set_step(MembraneStep.PREVIEW)
        st.rerun()

    render_protein_insertion(
        global_info,
        config.box_x, config.box_y, config.box_z, config.n_systems,
        config,
        temp_folder=temp_folder, protein_only=protein_only,
    )


def _validate(membrane_params) -> tuple[bool, list[str]]:
    errors = []
    entries = membrane_params["entries"]
    abs_vals = membrane_params["abs_lip_vals"]
    template_active = bool(membrane_params.get("template_active", False))

    if not entries:
        errors.append("No lipid entries defined.")
        return False, errors

    if isinstance(entries, dict):
        blocks = list(entries.items())
    else:
        blocks = [("membrane", entries)]

    for block_name, block_entries in blocks:
        if not block_entries:
            errors.append(f"Block '{block_name}': no lipid entries defined.")
            continue

        upper_vals = [e[1] for e in block_entries]
        lower_vals = [e[2] for e in block_entries]
        sum_upper = sum(upper_vals)
        sum_lower = sum(lower_vals)

        if abs_vals:
            if abs(sum_upper - round(sum_upper)) > 1e-6:
                errors.append(f"Block '{block_name}' upper leaflet: absolute values must sum to a whole number (got {sum_upper}).")
            if abs(sum_lower - round(sum_lower)) > 1e-6:
                errors.append(f"Block '{block_name}' lower leaflet: absolute values must sum to a whole number (got {sum_lower}).")
        elif not template_active:
            tolerance = 1e-6
            if abs(sum_upper - 1.0) > tolerance:
                errors.append(f"Block '{block_name}' upper leaflet: ratios must sum to 1.0 (got {sum_upper:.4f}).")
            if abs(sum_lower - 1.0) > tolerance:
                errors.append(f"Block '{block_name}' lower leaflet: ratios must sum to 1.0 (got {sum_lower:.4f}).")

    return len(errors) == 0, errors
