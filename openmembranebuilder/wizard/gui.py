import os
import tempfile
import numpy as np
import streamlit as st
from streamlit_molstar import st_molstar

_PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "NALA", "NARG", "NASN", "NASP", "NCYS", "NGLN", "NGLU",
    "NGLY", "NHIS", "NILE", "NLEU", "NLYS", "NMET", "NPHE",
    "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL",
    "CALA", "CARG", "CASN", "CASP", "CCYS", "CGLN", "CGLU",
    "CGLY", "CHIS", "CILE", "CLEU", "CLYS", "CMET", "CPHE",
    "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL",
}


def _prepare_display_pdb(pdb_path: str, box_nm):
    """
    Return a temporary PDB only when a periodic protein split can be repaired.

    The original simulation PDB is never modified.
    """
    if not pdb_path or box_nm is None or not os.path.exists(pdb_path):
        return pdb_path

    try:
        import MDAnalysis as mda

        universe = mda.Universe(pdb_path)

        protein = universe.select_atoms(
            "resname " + " ".join(sorted(_PROTEIN_RESNAMES))
        )

        if len(protein) == 0:
            return pdb_path

        box_lengths = np.asarray(box_nm, dtype=float) * 10.0
        universe.dimensions = np.asarray(
            [
                box_lengths[0],
                box_lengths[1],
                box_lengths[2],
                90.0,
                90.0,
                90.0,
            ],
            dtype=np.float32,
        )

        protein_positions = protein.positions.copy()
        original_span = protein_positions.max(axis=0) - protein_positions.min(axis=0)

        # If the protein is not split across a boundary, use the original PDB.
        if not np.any(original_span > 0.5 * box_lengths):
            return pdb_path

        repaired_positions = _repair_protein_periodic_split(
            protein_positions,
            box_lengths,
        )

        if repaired_positions is None:
            return pdb_path

        # Apply the repaired coordinates only to the protein atoms.
        protein.positions = repaired_positions

        session_dir = os.path.join(
            tempfile.gettempdir(),
            "openmembranebuilder-display",
        )
        os.makedirs(session_dir, exist_ok=True)

        display_path = os.path.join(
            session_dir,
            f"{os.path.basename(pdb_path)}-display-pbc-repaired.pdb",
        )

        with mda.Writer(display_path) as writer:
            writer.write(universe.atoms)

        return display_path

    except Exception as exc:
        st.warning(f"Could not repair PBC display for {pdb_path}: {exc}")
        return pdb_path

def _repair_protein_periodic_split(protein_positions, box_lengths):
    """
    Repair a single protein split across a periodic boundary.

    Returns repaired positions or None if no reliable repair is found.
    """
    repaired = protein_positions.copy()
    original_span = repaired.max(axis=0) - repaired.min(axis=0)

    for axis in range(3):
        span = original_span[axis]
        box_length = box_lengths[axis]

        # No likely PBC split in this direction.
        if span <= 0.5 * box_length:
            continue

        values = np.sort(repaired[:, axis])
        gaps = np.diff(values)

        if len(gaps) == 0:
            continue

        largest_gap_index = int(np.argmax(gaps))
        largest_gap = gaps[largest_gap_index]

        # The largest coordinate gap must plausibly be the periodic gap.
        if largest_gap <= 0.5 * box_length:
            continue

        split_value = values[largest_gap_index]

        # Atoms above the gap are shifted by -box length.
        candidate = repaired.copy()
        upper = candidate[:, axis] > split_value
        candidate[upper, axis] -= box_length

        candidate_span = candidate.max(axis=0) - candidate.min(axis=0)

        if candidate_span[axis] < span:
            repaired = candidate
            original_span = candidate_span

    if np.any(original_span > 0.5 * box_lengths):
        return None

    return repaired


def render_wizard_results(output_pdbs: list, zip_path: str, box_nm=None) -> None:
    st.title("Build Complete")

    if zip_path and os.path.exists(zip_path):
        with open(zip_path, "rb") as f:
            st.download_button(
                label=f"Download {os.path.basename(zip_path)}",
                data=f.read(),
                file_name=os.path.basename(zip_path),
                mime="application/zip",
            )

    st.divider()
    if not output_pdbs:
        st.warning("No output structures found.")
        return

    for i, (label, pdb_path) in enumerate(output_pdbs[:2]):
        display_path = _prepare_display_pdb(pdb_path, box_nm)
        with st.expander(f"System {i+1} — {label}", expanded=True):
            st_molstar(display_path, height=400, key=f"wizard-molstar-{label}")

    if len(output_pdbs) > 2:
        st.info(f"{len(output_pdbs) - 2} additional systems built not shown. Download the ZIP to access all.")
