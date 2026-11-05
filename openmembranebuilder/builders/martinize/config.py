from dataclasses import dataclass, field, asdict, fields
from typing import List, Dict, Any, Optional


@dataclass
class ProteinConfig:
    """Configuration settings for a single protein object."""
    name: str = ""  # chosen by the user; the files handed to the membrane builder are <name>.pdb and <name>.itp
    protein_input_mode: str = "upload_cg"  # "upload_cg" or "martinize"
    base_folder: str = ""  # work folder inside the session temp folder (set by the page)

    # Already coarse-grained protein uploads
    cg_pdb_path: str = ""
    cg_itp_path: str = ""
    cg_martini_version: str = "unknown"  # "unknown" | "Martini 2" | "Martini 3" (only for uploaded CG files)

    # Atomistic protein options
    atomistic_source: str = "Upload PDB"  # "Upload PDB" or "Fetch from PDB database"
    atomistic_pdb_path: str = ""
    fetch_pdb_id: str = ""

    # Preprocessing flags
    clean_structure: bool = True
    selected_chains: List[str] = field(default_factory=list)
    chain_ranges: Dict[str, List[int]] = field(default_factory=dict)

    # Mutation options
    do_mutation: bool = False
    mutation_chain: str = ""
    mutation_residue_idx: Optional[int] = None  # resid
    mutation_residue_name: str = ""
    mutation_target: str = ""  # target 3-letter amino acid code

    # Martinize2 options
    forcefield: str = "martini3001"  # "martini22", "martini3001", "martini3001-idp", etc.
    secondary_structure_mode: str = "MDTraj DSSP"  # "MDTraj DSSP", "All coil", "None", "Custom/Precalculated"
    custom_ss_string: str = ""
    use_elastic_network: bool = True
    elastic_lower: float = 0.5
    elastic_upper: float = 0.9
    elastic_force: float = 700.0
    position_restraints: str = "none"  # "none", "backbone", "all"
    cysteine_bridges: str = "auto"  # "auto", "detect", "none", etc.
    extra_flags: str = ""

    # Results, written by MartinizeRunner. These are what gets handed to the membrane
    # builder (single source of truth, replaces the ad-hoc session_state["protein_files"]).
    output_cg_pdb_path: str = ""
    output_itp_path: str = ""


def martini_ff_of(pc: "ProteinConfig") -> Optional[str]:
    """Membrane-builder force field ('martini_v2.2' / 'martini_v3') this protein belongs to,
    or None if unknown (uploaded CG files without a stated version)."""
    if pc.protein_input_mode == "martinize":
        ff = (pc.forcefield or "").lower()
        if ff.startswith("martini3"):
            return "martini_v3"
        if ff.startswith("martini2"):      # martini22, martini22p
            return "martini_v2.2"
        return None
    v = (pc.cg_martini_version or "").strip()
    if v.endswith("3"):
        return "martini_v3"
    if v.endswith("2"):
        return "martini_v2.2"
    return None


_PROTEIN_FIELDS = {f.name for f in fields(ProteinConfig)}


@dataclass
class Config:
    """Top-level configuration holding system-wide settings and protein entries."""
    build_mode: str = "membrane_only"  # "membrane_only" or "protein_membrane"
    n_proteins: int = 0
    proteins: List[ProteinConfig] = field(default_factory=list)

    def __post_init__(self):
        """Converts dictionary items in self.proteins into ProteinConfig instances."""
        converted_proteins: List[ProteinConfig] = []

        for item in self.proteins:
            if isinstance(item, dict):
                # ignore keys of older config versions (e.g. the removed copy_number)
                converted_proteins.append(ProteinConfig(**{k: v for k, v in item.items() if k in _PROTEIN_FIELDS}))
            elif isinstance(item, ProteinConfig):
                converted_proteins.append(item)
            else:
                raise TypeError(f"Unexpected type for protein entry: {type(item)}")

        self.proteins = converted_proteins
        self.n_proteins = len(self.proteins)

    def to_dict(self) -> Dict[str, Any]:
        """Recursively converts the dataclass hierarchy into a dictionary."""
        return asdict(self)

    def protein_file_names(self) -> List[Optional[str]]:
        """Name each protein gets in the membrane builder (file name of its CG pdb),
        aligned with `self.proteins`; None for proteins without results."""
        import os
        names: List[Optional[str]] = []
        used = set()
        for i, pc in enumerate(self.proteins):
            if not pc.output_cg_pdb_path or not pc.output_itp_path:
                names.append(None)
                continue
            name = os.path.basename(pc.output_cg_pdb_path)
            if name in used:  # same file name twice -> keep both, unique names
                name = f"p{i}_{name}"
            used.add(name)
            names.append(name)
        return names

    def protein_files(self) -> Dict[str, Dict[str, str]]:
        """{pdb_name: {pdb_path, itp_path, force_field}} for the membrane builder, built from
        the CURRENT proteins only (no stale entries after delete / re-run / back).
        force_field is 'martini_v2.2', 'martini_v3' or '' (unknown)."""
        files: Dict[str, Dict[str, str]] = {}
        for name, pc in zip(self.protein_file_names(), self.proteins):
            if name:
                files[name] = {
                    "pdb_path": pc.output_cg_pdb_path,
                    "itp_path": pc.output_itp_path,
                    "force_field": martini_ff_of(pc) or "",
                }
        return files
