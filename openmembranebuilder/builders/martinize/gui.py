import os
from pathlib import Path
from typing import Dict, Any, List
import streamlit as st

from openmembranebuilder.builders.martinize.config import Config, ProteinConfig


# ─── Page 1: Protein Input ────────────────────────────────────────────────────

def render_input(config: Config, temp_folder: str) -> None:
    """Render the protein input page.

    Shows one panel per protein in config.proteins.  Each panel has its own
    header ("Protein N"), a delete button, mode selection, file upload / RCSB
    fetch and a clean-structure checkbox.

    All widget keys are namespaced with ``_p{i}`` so panels don't collide.

    Sets the following session_state signals:
    - ``_martinize_input_action``  : "next" | "back"
    - ``_martinize_add_protein``   : True  (add a blank protein)
    - ``_martinize_delete_protein``: int   (index of protein to delete)
    - ``_martinize_fetch_requested``: int  (index of protein to fetch PDB for)
    """
    st.title("🧬 Protein Structure Input")
    st.caption("Provide the protein structures to be coarse-grained with martinize2.")

    temp_dir = Path(temp_folder) / "protein_input"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # ── One panel per protein ─────────────────────────────────────────────────
    for i, pconfig in enumerate(config.proteins):
        _render_protein_panel(i, pconfig, config.n_proteins, temp_dir)

    # ── Add a protein ─────────────────────────────────────────────────────────
    st.divider()
    if st.button("➕ Add a protein", use_container_width=False):
        st.session_state["_martinize_add_protein"] = True
        st.rerun()

    # ── Navigation ────────────────────────────────────────────────────────────
    with st.container(key="nav_martinize_input"):
        st.divider()
        left, _, right = st.columns([1, 6, 1])
        if left.button("← Back", use_container_width=True):
            st.session_state["_martinize_input_action"] = "back"
            st.rerun()
        if right.button("Next →", use_container_width=True, type="primary"):
            st.session_state["_martinize_input_action"] = "next"
            st.rerun()


def _render_protein_panel(
    i: int,
    pconfig: ProteinConfig,
    n_proteins: int,
    temp_dir: Path,
) -> None:
    """Render the input panel for the i-th protein (0-based index)."""
    with st.container(border=True):
        # ── Header row ────────────────────────────────────────────────────────
        hdr_left, hdr_right = st.columns([4, 1])
        hdr_left.subheader(f"Protein {i + 1}")
        if hdr_right.button(
            "🗑 Delete", key=f"del_protein_{i}",
            use_container_width=True,
            help="Remove this protein from the workflow",
        ):
            st.session_state["_martinize_delete_protein"] = i
            st.rerun()

        # ── Name: the files become <name>.pdb and <name>.itp ──────────────────
        name_key = f"protein_name_p{i}"
        st.session_state.setdefault(name_key, pconfig.name)
        st.text_input(
            "Protein name",
            key=name_key,
            placeholder="e.g. 1UBQ",
            help="Required. The structure and topology files are named <name>.pdb and <name>.itp "
                 "(letters, digits, _ and -; each protein needs its own name).",
        )

        # ── Mode radio ────────────────────────────────────────────────────────
        mode_key = f"protein_input_mode_p{i}"
        protein_input_mode = st.radio(
            "How do you want to provide the protein structure?",
            ["upload_cg", "martinize"],
            format_func=lambda x: {
                "upload_cg": "Upload already coarse-grained PDB/GRO + topology ITP",
                "martinize": "Start from atomistic PDB (Download from RCSB or Upload)",
            }[x],
            index=["upload_cg", "martinize"].index(
                st.session_state.get(mode_key, pconfig.protein_input_mode)
            ),
            key=mode_key,
        )

        st.divider()

        if protein_input_mode == "upload_cg":
            _render_upload_cg_panel(i, pconfig, temp_dir)
        else:
            _render_atomistic_panel(i, pconfig, temp_dir)



def _render_upload_cg_panel(i: int, pconfig: ProteinConfig, temp_dir: Path) -> None:
    st.subheader("Upload Coarse-Grained Protein Files")
    col1, col2 = st.columns(2)

    with col1:
        uploaded_cg_pdb = st.file_uploader(
            "Upload CG Structure (.pdb, .gro)", type=["pdb", "gro"],
            key=f"cg_pdb_uploader_p{i}",
        )
        if uploaded_cg_pdb:
            dest = str(temp_dir / f"p{i}_{uploaded_cg_pdb.name}")
            with open(dest, "wb") as f:
                f.write(uploaded_cg_pdb.getvalue())
            st.session_state[f"cg_pdb_path_p{i}"] = dest
            st.success(f"Uploaded: {uploaded_cg_pdb.name}")
        elif f"cg_pdb_path_p{i}" not in st.session_state:
            st.session_state[f"cg_pdb_path_p{i}"] = pconfig.cg_pdb_path

    with col2:
        uploaded_cg_itp = st.file_uploader(
            "Upload Protein Topology (.itp)", type=["itp"],
            key=f"cg_itp_uploader_p{i}",
        )
        if uploaded_cg_itp:
            dest = str(temp_dir / f"p{i}_{uploaded_cg_itp.name}")
            with open(dest, "wb") as f:
                f.write(uploaded_cg_itp.getvalue())
            st.session_state[f"cg_itp_path_p{i}"] = dest
            st.success(f"Uploaded: {uploaded_cg_itp.name}")
        elif f"cg_itp_path_p{i}" not in st.session_state:
            st.session_state[f"cg_itp_path_p{i}"] = pconfig.cg_itp_path

    version_options = ["unknown", "Martini 2", "Martini 3"]
    version_key = f"cg_version_p{i}"
    st.session_state.setdefault(version_key, pconfig.cg_martini_version if pconfig.cg_martini_version in version_options else "unknown")
    st.selectbox(
        "Martini version of these files",
        version_options,
        key=version_key,
        help="Used to warn you if the protein does not match the force field of the membrane. "
             "Martini 2 and Martini 3 models must not be mixed.",
    )


def _render_atomistic_panel(i: int, pconfig: ProteinConfig, temp_dir: Path) -> None:
    st.subheader("Atomistic Structure Input")

    source_key = f"atomistic_source_p{i}"
    atomistic_source = st.radio(
        "Select atomistic structure source:",
        ["Upload PDB", "Fetch from PDB database"],
        index=["Upload PDB", "Fetch from PDB database"].index(
            st.session_state.get(source_key, pconfig.atomistic_source)
        ),
        key=source_key,
    )

    if atomistic_source == "Upload PDB":
        uploaded_at = st.file_uploader(
            "Upload Atomistic PDB File", type=["pdb"],
            key=f"atomistic_pdb_uploader_p{i}",
        )
        if uploaded_at:
            dest = str(temp_dir / f"p{i}_{uploaded_at.name}")
            with open(dest, "wb") as f:
                f.write(uploaded_at.getvalue())
            prev = st.session_state.get(f"atomistic_path_p{i}", "")
            if prev != dest:
                _reset_downstream_state(i)
            st.session_state[f"atomistic_path_p{i}"] = dest
            st.success(f"Uploaded: {uploaded_at.name}")
    else:
        col1, col2 = st.columns([3, 1])
        with col1:
            st.text_input(
                "Enter 4-letter RCSB PDB ID:",
                value=st.session_state.get(f"fetch_pdb_id_p{i}", pconfig.fetch_pdb_id),
                placeholder="e.g. 1UBQ",
                key=f"fetch_pdb_id_p{i}",
            )
        with col2:
            st.write(" ")
            st.write(" ")
            if st.button("Fetch Structure", type="primary", key=f"fetch_btn_p{i}"):
                pdb_id = st.session_state.get(f"fetch_pdb_id_p{i}", "").strip()
                if len(pdb_id) == 4:
                    st.session_state["_martinize_fetch_requested"] = i
                    st.rerun()
                else:
                    st.error("PDB ID must be exactly 4 characters.")

        if st.session_state.get(f"fetched_pdb_success_p{i}"):
            st.success(st.session_state[f"fetched_pdb_success_p{i}"])

    if st.session_state.get(f"atomistic_path_p{i}"):
        st.write(f"**Loaded structure:** `{st.session_state[f'atomistic_path_p{i}']}`")

    # Clean-up checkbox
    st.divider()
    st.subheader("Structure Clean-up")
    st.checkbox(
        "Clean structure with MDAnalysis (removes water, ligands, ions – keeps protein only)",
        value=st.session_state.get(f"clean_structure_p{i}", pconfig.clean_structure),
        key=f"clean_structure_p{i}",
    )


# ─── Page 2: Configure & Martinize ────────────────────────────────────────────

def render_configure(
    pconfig: ProteinConfig,
    protein_idx: int,
    chains: List[str],
    summaries: Dict[str, Any],
) -> None:
    """Render the configuration page for a single protein.

    Pure presentation only.  Sets ``_martinize_configure_action`` to
    "back", "run", or "dssp_preview" on button presses.

    All widget keys are namespaced with ``_p{protein_idx}`` so state is
    independent for each protein.
    """
    p = protein_idx  # shorthand for key namespacing

    st.title("⚙️ Configure & Martinize")
    st.caption(f"Protein {protein_idx + 1} — Select chains, configure mutations, and set martinize2 parameters.")

    # ── Chain selection ───────────────────────────────────────────────────────
    st.subheader("Chain Selection")

    new_selected: List[str] = []
    new_ranges: Dict[str, List[int]] = {}

    if chains:
        prev_selected: List[str] = st.session_state.get(
            f"selected_chains_p{p}", pconfig.selected_chains or chains
        )
        prev_selected = [c for c in prev_selected if c in chains]
        chain_ranges: Dict[str, List[int]] = st.session_state.get(
            f"chain_ranges_p{p}", dict(pconfig.chain_ranges)
        )

        hdr = st.columns([0.5, 1, 1.5, 1.5, 1.5, 1.5])
        hdr[0].markdown("**Select**")
        hdr[1].markdown("**Chain ID**")
        hdr[2].markdown("**# Residues**")
        hdr[3].markdown("**First ResID**")
        hdr[4].markdown("**Last ResID**")
        hdr[5].markdown("**Residue Range**")

        for c in chains:
            s = summaries.get(c, {})
            row = st.columns([0.5, 1, 1.5, 1.5, 1.5, 1.5])
            is_selected = row[0].checkbox(
                "sel", value=(c in prev_selected),
                key=f"chain_sel_{c}_p{p}", label_visibility="collapsed",
            )
            row[1].write(c)
            if s:
                row[2].write(str(s.get("n_residues", "—")))
                row[3].write(str(s.get("first_resid", "—")))
                row[4].write(str(s.get("last_resid", "—")))
            else:
                for col in row[2:5]:
                    col.write("—")

            if is_selected:
                new_selected.append(c)
                if s:
                    first_resid = s.get("first_resid", 1)
                    last_resid  = s.get("last_resid", 9999)
                    prev_range  = chain_ranges.get(c, [first_resid, last_resid])
                    prev_start  = prev_range[0] if prev_range and len(prev_range) == 2 else first_resid
                    prev_end    = prev_range[1] if prev_range and len(prev_range) == 2 else last_resid
                    rc = row[5].columns(2)
                    range_start = rc[0].number_input(
                        "From", min_value=first_resid, max_value=last_resid,
                        value=max(first_resid, min(last_resid, prev_start)),
                        key=f"chain_range_start_{c}_p{p}", label_visibility="collapsed",
                    )
                    range_end = rc[1].number_input(
                        "To", min_value=first_resid, max_value=last_resid,
                        value=max(first_resid, min(last_resid, prev_end)),
                        key=f"chain_range_end_{c}_p{p}", label_visibility="collapsed",
                    )
                    new_ranges[c] = [int(range_start), int(range_end)]
                else:
                    new_ranges[c] = []

        if not new_selected:
            st.warning("⚠️ Please select at least one chain.")

        st.session_state[f"selected_chains_p{p}"] = new_selected
        st.session_state[f"chain_ranges_p{p}"]    = new_ranges

    else:
        st.warning("No protein chains could be identified in the structure.")

    # ── Mutation ──────────────────────────────────────────────────────────────
    st.divider()
    st.subheader("Mutation")
    do_mutation = st.checkbox(
        "Introduce residue mutation",
        value=st.session_state.get(f"do_mutation_p{p}", pconfig.do_mutation),
        key=f"do_mutation_p{p}",
    )

    if do_mutation and new_selected:
        mutation_chain       = st.session_state.get(f"mutation_chain_p{p}",       pconfig.mutation_chain)
        mutation_residue_idx = st.session_state.get(f"mutation_residue_idx_p{p}", pconfig.mutation_residue_idx)
        mutation_target      = st.session_state.get(f"mutation_target_p{p}",      pconfig.mutation_target)

        if mutation_chain not in new_selected:
            mutation_chain = new_selected[0]

        hdr = st.columns([1.2, 1.2, 1.2, 1.2])
        hdr[0].markdown("**Chain ID**")
        hdr[1].markdown("**Amino acid**")
        hdr[2].markdown("**ResID**")
        hdr[3].markdown("**Mutant**")

        sel = st.columns([1.2, 1.2, 1.2, 1.2])

        mutation_chain = sel[0].selectbox(
            "Chain", options=new_selected,
            index=new_selected.index(mutation_chain),
            key=f"mut_chain_sel_p{p}", label_visibility="collapsed",
        )

        all_residues = st.session_state.get(f"residues_{mutation_chain}_p{p}", [])
        chosen_range = new_ranges.get(mutation_chain, [])
        if chosen_range and len(chosen_range) == 2:
            range_start_val, range_end_val = chosen_range
            range_residues = [r for r in all_residues if range_start_val <= r["resid"] <= range_end_val]
        else:
            range_residues = all_residues

        if range_residues:
            unique_names = sorted({r["resname"] for r in range_residues})
            aa_filter_opts = ["All"] + unique_names
            prev_aa_filter = st.session_state.get(f"mut_aa_filter_p{p}", "All")
            if prev_aa_filter not in aa_filter_opts:
                prev_aa_filter = "All"
            aa_filter = sel[1].selectbox(
                "AA filter", options=aa_filter_opts,
                index=aa_filter_opts.index(prev_aa_filter),
                key=f"mut_aa_filter_p{p}", label_visibility="collapsed",
            )

            visible_residues = range_residues if aa_filter == "All" else [
                r for r in range_residues if r["resname"] == aa_filter
            ]
            resid_list = [r["resid"] for r in visible_residues]

            if resid_list:
                if mutation_residue_idx not in resid_list:
                    mutation_residue_idx = resid_list[0]
                sel[2].selectbox(
                    "RESID", options=resid_list,
                    index=resid_list.index(mutation_residue_idx),
                    key=f"mut_resid_sel_p{p}", label_visibility="collapsed",
                )
                standard_aas = [
                    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                    "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                ]
                sel[3].selectbox(
                    "Mutant", options=standard_aas,
                    index=standard_aas.index(mutation_target) if mutation_target in standard_aas else 0,
                    key=f"mut_target_sel_p{p}", label_visibility="collapsed",
                )
            else:
                st.warning(f"No residues of type **{aa_filter}** found in the selected range.")
        else:
            st.warning("No residues found in the selected range for this chain.")
    elif do_mutation and not new_selected:
        st.warning("⚠️ Select at least one chain above before configuring a mutation.")

    # ── Martinize2 Options ────────────────────────────────────────────────────
    st.divider()
    st.subheader("Martinize Coarse-Graining Options")

    forcefield_opts = ["martini22", "martini3001", "martini3001-idp"]
    st.selectbox(
        "Force field (-ff):",
        forcefield_opts,
        index=forcefield_opts.index(st.session_state.get(f"forcefield_p{p}", pconfig.forcefield)),
        key=f"forcefield_p{p}",
    )

    ss_opts = ["MDTraj DSSP", "All coil", "None", "Custom/Precalculated"]
    st.selectbox(
        "Secondary structure mode:",
        ss_opts,
        index=ss_opts.index(st.session_state.get(f"secondary_structure_mode_p{p}", pconfig.secondary_structure_mode)),
        key=f"secondary_structure_mode_p{p}",
    )

    if st.session_state.get(f"secondary_structure_mode_p{p}") == "Custom/Precalculated":
        st.text_input(
            "Precalculated secondary structure string (-ss):",
            value=st.session_state.get(f"custom_ss_string_p{p}", pconfig.custom_ss_string),
            placeholder="e.g. HHHHCCCCCCEEEEEE",
            key=f"custom_ss_string_p{p}",
        )

    if st.session_state.get(f"secondary_structure_mode_p{p}") == "MDTraj DSSP":
        if st.button("Preview calculated DSSP sequence", key=f"dssp_btn_p{p}"):
            st.session_state["_martinize_configure_action"] = "dssp_preview"
            st.rerun()
        if f"dssp_preview_result_p{p}" in st.session_state:
            st.text_area(
                "DSSP String Preview:", value=st.session_state[f"dssp_preview_result_p{p}"],
                height=100, disabled=True, key=f"dssp_area_p{p}",
            )
            if f"dssp_preview_error_p{p}" in st.session_state:
                st.error(st.session_state.pop(f"dssp_preview_error_p{p}"))

    st.checkbox(
        "Enable elastic network (-elastic)",
        value=st.session_state.get(f"use_elastic_network_p{p}", pconfig.use_elastic_network),
        key=f"use_elastic_network_p{p}",
    )
    if st.session_state.get(f"use_elastic_network_p{p}", pconfig.use_elastic_network):
        c1, c2, c3 = st.columns(3)
        c1.number_input(
            "Lower cutoff (-el nm)", min_value=0.0, max_value=5.0, step=0.1,
            value=st.session_state.get(f"elastic_lower_p{p}", pconfig.elastic_lower),
            key=f"elastic_lower_p{p}",
        )
        c2.number_input(
            "Upper cutoff (-eu nm)", min_value=0.0, max_value=5.0, step=0.1,
            value=st.session_state.get(f"elastic_upper_p{p}", pconfig.elastic_upper),
            key=f"elastic_upper_p{p}",
        )
        c3.number_input(
            "Force constant (-ef)", min_value=0.0, max_value=10000.0, step=100.0,
            value=st.session_state.get(f"elastic_force_p{p}", pconfig.elastic_force),
            key=f"elastic_force_p{p}",
        )

    restraints_opts = ["none", "backbone", "all"]
    st.selectbox(
        "Position restraints (-p):",
        restraints_opts,
        index=restraints_opts.index(st.session_state.get(f"position_restraints_p{p}", pconfig.position_restraints)),
        key=f"position_restraints_p{p}",
    )

    cys_opts = ["auto", "detect", "none"]
    st.selectbox(
        "Cysteine bridges (-cys):",
        cys_opts,
        index=cys_opts.index(st.session_state.get(f"cysteine_bridges_p{p}", pconfig.cysteine_bridges)),
        key=f"cysteine_bridges_p{p}",
    )

    st.text_input(
        "Extra martinize2 options (space-separated):",
        value=st.session_state.get(f"extra_flags_p{p}", pconfig.extra_flags),
        placeholder="e.g. -maxwarn 5 -nt",
        key=f"extra_flags_p{p}",
    )

    # ── Navigation ────────────────────────────────────────────────────────────
    with st.container(key="nav_martinize_configure"):
        st.divider()
        left, mid, right = st.columns([1, 3, 1])
        if left.button("← Back", use_container_width=True):
            st.session_state["_martinize_configure_action"] = "back"
            st.rerun()
        if mid.button("🧬 Run Martinize Preprocessor", use_container_width=True, type="primary"):
            st.session_state["_martinize_configure_action"] = "run"
            st.rerun()


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _reset_downstream_state(protein_idx: int) -> None:
    """Clear all session state that depends on the loaded structure for protein i."""
    p = protein_idx
    for key in [
        f"cleaned_path_p{p}", f"fetched_pdb_success_p{p}",
        f"selected_chains_p{p}", f"chain_ranges_p{p}",
        f"do_mutation_p{p}", f"mutation_chain_p{p}",
        f"mutation_residue_idx_p{p}", f"mutation_residue_name_p{p}",
        f"mutation_target_p{p}", f"dssp_preview_result_p{p}",
        f"active_struct_p{p}", f"chains_p{p}", f"chain_summaries_p{p}",
    ]:
        st.session_state.pop(key, None)
    for key in list(st.session_state.keys()):
        if key.startswith(f"residues_") and key.endswith(f"_p{p}"):
            del st.session_state[key]
