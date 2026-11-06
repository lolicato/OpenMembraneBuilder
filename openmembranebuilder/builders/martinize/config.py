"""The settings of the protein preparation.

There is one `ProteinConfig` per protein and one `Config` that holds them all. When the
proteins are ready, `Config.protein_files` hands their names and files to the membrane
builder.
"""

from dataclasses import dataclass, field, asdict, fields
from typing import List, Dict, Any, Optional


@dataclass
class ProteinConfig:
    """Everything about one protein: its name, where it comes from and how it is prepared.

    Attributes:
        name: The name you choose (e.g. `1UBQ`). The files handed to the membrane builder
            are `<name>.pdb` and `<name>.itp`, and the topology lists the protein under
            this name. Letters, digits, `_` and `-`; each protein needs its own name.
        protein_input_mode: `"upload_cg"` (already coarse-grained files) or `"martinize"`
            (start from an atomistic structure).
        base_folder: Work folder of this protein inside the session folder.
        cg_pdb_path: Uploaded coarse-grained structure (`.pdb` or `.gro`).
        cg_itp_path: Uploaded coarse-grained topology (`.itp`).
        cg_martini_version: Martini generation of uploaded files: `"unknown"`,
            `"Martini 2"` or `"Martini 3"` (used to warn about mixing generations).
        atomistic_source: `"Upload PDB"` or `"Fetch from PDB database"`.
        atomistic_pdb_path: The atomistic structure in use.
        fetch_pdb_id: Four-character RCSB PDB code to download, e.g. `1UBQ`.
        clean_structure: Keep only the protein (remove water, ligands, ions).
        selected_chains: The chains to keep.
        chain_ranges: Residue range to keep per chain: `{"A": [first, last]}`.
        do_mutation: Mutate one residue before coarse-graining (needs PyMOL).
        mutation_chain: Chain of the residue to mutate.
        mutation_residue_idx: Residue number to mutate.
        mutation_residue_name: Current amino acid of that residue.
        mutation_target: New amino acid (three-letter code, e.g. `ALA`).
        forcefield: martinize2 force field: `"martini3001"`, `"martini3001-idp"` or
            `"martini22"`.
        secondary_structure_mode: Where the secondary structure comes from:
            `"MDTraj DSSP"` (calculated), `"All coil"`, `"None"` or
            `"Custom/Precalculated"`.
        custom_ss_string: Your own secondary structure, one letter per residue.
        use_elastic_network: Add an elastic network that keeps the fold.
        elastic_lower: Elastic network lower cut-off, in nm.
        elastic_upper: Elastic network upper cut-off, in nm.
        elastic_force: Elastic network force constant, in kJ mol⁻¹ nm⁻².
        position_restraints: Position restraints written by martinize2: `"none"`,
            `"backbone"` or `"all"`.
        cysteine_bridges: Disulfide bridges: `"auto"`, `"detect"` or `"none"`.
        extra_flags: Further martinize2 options, separated by spaces.
        output_cg_pdb_path: The finished coarse-grained structure (`<name>.pdb`).
        output_itp_path: The finished topology (`<name>.itp`).
    """
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
    """The membrane force field a protein belongs to.

    Args:
        pc: The protein's settings.

    Returns:
        `"martini_v3"`, `"martini_v2.2"`, or `None` when unknown (uploaded files without a
        stated generation).
    """
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
    """The settings of all proteins of the build.

    Attributes:
        build_mode: `"protein_membrane"` when proteins are prepared (always so in the
            wizard), `"membrane_only"` otherwise.
        n_proteins: Number of proteins.
        proteins: The settings of every protein.
    """
    build_mode: str = "membrane_only"  # "membrane_only" or "protein_membrane"
    n_proteins: int = 0
    proteins: List[ProteinConfig] = field(default_factory=list)

    def __post_init__(self):
        """Turn proteins read from a JSON file into `ProteinConfig` objects.

        Settings that older versions of the program stored but no longer use are ignored.

        Raises:
            TypeError: A protein entry is neither a dictionary nor a `ProteinConfig`.
        """
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
        """All settings as plain dictionaries and lists (for JSON)."""
        return asdict(self)

    def protein_file_names(self) -> List[Optional[str]]:
        """The file name each finished protein gets in the membrane builder.

        Returns:
            One name per protein (e.g. `"1UBQ.pdb"`), in the order of `proteins`; `None` for
            proteins that are not finished. Should two files have the same name, the second is
            renamed so both are kept.
        """
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
        """The finished proteins, as the membrane builder expects them.

        Only the current proteins are listed, so a deleted or re-run protein never leaves an
        old entry behind.

        Returns:
            `{"1UBQ.pdb": {"pdb_path", "itp_path", "force_field"}}`; `force_field` is
            `"martini_v3"`, `"martini_v2.2"` or empty when unknown.
        """
        files: Dict[str, Dict[str, str]] = {}
        for name, pc in zip(self.protein_file_names(), self.proteins):
            if name:
                files[name] = {
                    "pdb_path": pc.output_cg_pdb_path,
                    "itp_path": pc.output_itp_path,
                    "force_field": martini_ff_of(pc) or "",
                }
        return files
