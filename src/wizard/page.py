import os
import json
import random
import shutil
from datetime import datetime
from collections import OrderedDict
from dataclasses import asdict, is_dataclass

import streamlit as st

from src.wizard.state import get_state, WizardStep, PATHWAYS
from src.wizard.navigation import go_back, go_forward
from src.core.global_info import GlobalInfo


# ─── Temp folder handling ─────────────────────────────────────────────────────

def _ensure_temp_folder(state) -> str:
    if not state.temp_folder:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        rand = random.randint(100000, 999999)
        state.temp_folder = f"./temp-uploads_{ts}-{rand}"
    return state.temp_folder


def _cleanup_temp_folder(state) -> None:
    if state.temp_folder and os.path.exists(state.temp_folder):
        try:
            shutil.rmtree(state.temp_folder)
        except Exception as e:
            print(f"Could not remove temp folder {state.temp_folder}: {e}")
    state.temp_folder = ""


# ─── Config export ────────────────────────────────────────────────────────────

def _json_fallback(obj):
    if hasattr(obj, "__dict__"):
        return vars(obj)
    if hasattr(obj, "__iter__"):
        return list(obj)
    return str(obj)


# Keys that are transient / machine specific.
# NOTE: "molecule_builder_index" is intentionally NOT stripped any more:
# without it a headless run cannot rebuild the custom lipids.
_STRIP_KEYS = {
    "lipid_entries_absolute",
    "lipid_entries_relative",
    "membrane_string",
    "entries_current",
    "base_folder",
    "temp_folder",
    "external_lipid_definitions",
    "validation_output_path",
    "metadata_output_path",
    "metadata",
}

_ORDER = ["martinize", "membrane"]


def _save_full_config(config_map: dict, output_folder: str) -> str:
    combined = OrderedDict()
    for key in _ORDER:
        if key not in config_map:
            continue
        cfg = config_map[key]
        data = asdict(cfg) if is_dataclass(cfg) else vars(cfg).copy()
        for k in _STRIP_KEYS:
            data.pop(k, None)
        combined[key] = data

    os.makedirs(output_folder, exist_ok=True)
    path = os.path.join(output_folder, "config.json")
    with open(path, "w") as f:
        json.dump(combined, f, indent=4, default=_json_fallback)
    return path



# ─── Main render ──────────────────────────────────────────────────────────────

def render_wizard():
    state = get_state()

    if state.locked:
        with st.sidebar:
            if st.button("↩ Change workflow", use_container_width=True):
                st.session_state.clear()
                st.rerun()

    if state.step == WizardStep.SELECT_PATHWAY:
        _select_pathway(state)

    elif state.step == WizardStep.MARTINIZE:
        from src.builders.martinize.page import main as martinize_main

        result, cfg = martinize_main(temp_folder=_ensure_temp_folder(state))
        if result == "done" and cfg:
            st.session_state.setdefault("completed_configs", {})["martinize"] = cfg
            go_forward(state)
            st.rerun()
        elif result == "back":
            go_back(state)
            st.rerun()

    elif state.step == WizardStep.MEMBRANE:
        _render_membrane_stage(state)

    elif state.step == WizardStep.RESULTS:
        from src.wizard.gui import render_wizard_results

        membrane_cfg = st.session_state.get("completed_configs", {}).get("membrane")
        box_nm = None
        if membrane_cfg is not None:
            box_nm = (
                float(membrane_cfg.box_x),
                float(membrane_cfg.box_y),
                float(membrane_cfg.box_z),
            )

        render_wizard_results(
            st.session_state.get("output_pdbs", []),
            st.session_state.get("zip_path"),
            box_nm=box_nm,
        )


def _render_membrane_stage(state) -> None:
    completed = st.session_state.setdefault("completed_configs", {})

    from src.builders.membrane.page import main as membrane_main


   

    result, cfg = membrane_main(
        GlobalInfo(),
        has_protein=state.has_protein,
        temp_folder=_ensure_temp_folder(state),
        protein_only=state.pathway.protein_only if state.pathway else False
    )

    if result == "done" and cfg:
        completed["membrane"] = cfg
        # output_pdbs / zip_path are already stored in session_state by the membrane page
        config_path = _save_full_config(completed, cfg.base_folder)
        st.session_state["final_config_path"] = config_path
        _cleanup_temp_folder(state)
        go_forward(state)
        st.rerun()

    elif result == "back":
        go_back(state)
        st.rerun()


def _select_pathway(state) -> None:
    st.title("Welcome")
    st.subheader("Choose your workflow")

    labels = [p.label for p in PATHWAYS]
    default_label = state.pathway.label if state.pathway else "Membrane only (no protein)"

    chosen_label = st.radio(
        "What do you want to build?",
        labels,
        index=labels.index(default_label),
    )
    selected = next(p for p in PATHWAYS if p.label == chosen_label)

    st.divider()


    if st.button("Start →", type="primary", use_container_width=True):
        # a fresh start must not inherit results of a previous run
        st.session_state.pop("completed_configs", None)
        st.session_state.pop("protein_files", None)
        state.pathway = selected
        state.locked = True
        go_forward(state)
        st.rerun()
