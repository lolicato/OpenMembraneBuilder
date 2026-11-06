"""The settings of the membrane builder.

All choices made on the membrane, protein-insertion and preview pages end up in one `Config`
object. It is what the runner builds from, and (without the session-only values) what is
saved as the `membrane` section of `config.json`.
"""

from dataclasses import dataclass, asdict, field
import json


@dataclass
class Config:
    """Everything the membrane builder needs to build the systems.

    Attributes:
        output_name: Optional name added to the output folder
            (`OpenMembraneBuilder-<date>_<time>-<output_name>`).
        selected_ff: Force field of the membrane, `"martini_v3"` or `"martini_v2.2"`.
        base_folder: The output folder of this build (set when the build starts).
        box_x: Box length along x, in nm.
        box_y: Box length along y, in nm (rectangular boxes only).
        box_z: Box height along z, in nm (not used by skewed hexagonal boxes).
        box_type: `"rectangular"`, `"hexagonal"` or `"skewed_hexagonal"`.
        solvation: Water and ions in COBY's notation, e.g.
            `"solv:W pos:NA neg:CL salt_molarity:0.15"`.
        salt_molarity: Salt concentration in mol/L (also contained in `solvation`).
        n_systems: Number of independent systems (replicas) to build.
        protein_files: The proteins to insert:
            `{"1UBQ.pdb": {"pdb_path": ..., "itp_path": ..., "force_field": ...}}`.
        template_active: The composition comes from a membrane template.
        template_path: The template file in use.
        membrane_string: The composition in COBY's notation (made automatically).
        abs_lip_vals: `True`: the leaflet values are numbers of lipids;
            `False`: they are fractions that add up to 1 in each leaflet.
        lipid_entries_relative: Composition kept while in fraction mode (session only).
        lipid_entries_absolute: Composition kept while in number mode (session only).
        entries: The composition. Either one list of rows
            `[lipid, upper, lower, area per lipid upper, area per lipid lower]`, or, with a
            template that has named blocks, `{block name: list of rows}`.
        entries_current: The composition being built at the moment (session only).
        z_method: How the protein height is given: `"Absolute z position"` (z measured
            from the box centre) or `"Height above Membrane"` (gap between membrane surface
            and the lowest protein bead).
        distance_to_mem: The gap for "Height above Membrane", in nm.
        randomize_pos: Place the protein at a random x/y position.
        randomize_pos_every: Draw a new random position for every system.
        randomize_rot: Give the protein a random orientation.
        randomize_rot_every: Draw a new random orientation for every system.
        more_copies: Insert several copies of one protein.
        copy_number: How many copies per system.
        copy_mode: `"together"`: the copies sit on a regular grid and share one
            orientation; `"independent"`: every copy has its own position and orientation.
        copy_params: The pose of each independent copy:
            `[{"cx", "cy", "cz", "rx", "ry", "rz"}, ...]`.
        randomize_rot_all_copies: Give every copy ("together" mode) its own random
            orientation.
        protein_params: The pose actually used for every protein in every system:
            `{"R0001": {"1UBQ.pdb": {"cx", "cy", "cz", "rx", "ry", "rz"}}, ...}`.
            Positions in nm from the box centre, rotations in degrees.
        insertion_params: All protein-insertion settings per protein, as entered on the
            page (used to restore the page).
        use_protein_distance: Two proteins: place them by the gap between them instead of
            by x/y positions.
        protein_distance: That gap, in nm.
        molecule_builder_index: Custom lipids made by a molecule builder, by lipid name.
        external_lipid_definitions: Extra lipid descriptions supplied from outside
            (session only).
        extern_config: Settings handed in by another tool (session only).
    """

    output_name: str = ""
    selected_ff: str = "martini_v3"
    base_folder: str = ""

    box_x: float = 10.0
    box_y: float = 10.0
    box_z: float = 10.0
    box_type: str = "rectangular"

    solvation: str = ""
    salt_molarity: float = 0.15

    n_systems: int = 1

    protein_files: dict = field(default_factory=dict)

    template_active: bool = False
    template_path: str = ""
    membrane_string: str = ""

    abs_lip_vals: bool = False
    lipid_entries_relative: list = None
    lipid_entries_absolute: list = None
    entries: list = None
    entries_current: list = None  # transient: set by the runner per template

    z_method: str = "Absolute z position"
    distance_to_mem: float = 2.0
    randomize_pos: bool = False
    randomize_pos_every: bool = False
    randomize_rot: bool = False
    randomize_rot_every: bool = False
    more_copies: bool = False
    copy_number: int = 2
    copy_mode: str = "together"  # "together": copies spread on a grid | "independent": one pose per copy
    copy_params: list = field(default_factory=list)  # independent copies: [{cx, cy, cz, rx, ry, rz}, ...]
    randomize_rot_all_copies: bool = False
    protein_params: dict = field(default_factory=dict, repr=False, compare=False)

    # previously ad-hoc attributes -> now real fields so asdict()/JSON keep them
    insertion_params: dict = field(default_factory=dict, repr=False, compare=False)
    use_protein_distance: bool = False
    protein_distance: float = 1.0

    molecule_builder_index: dict = field(default_factory=dict, repr=False, compare=False)
    external_lipid_definitions: dict = field(default_factory=dict, repr=False, compare=False)
    extern_config: dict = field(default_factory=dict, repr=False, compare=False)  # transient

    def __post_init__(self):
        """Fill in empty lists and give every system a starting pose (all zeros) for the proteins.
        """
        if self.lipid_entries_relative is None:
            self.lipid_entries_relative = []
        if self.lipid_entries_absolute is None:
            self.lipid_entries_absolute = []
        if self.entries is None:
            self.entries = []
        if self.entries_current is None:
            self.entries_current = []
        if not self.protein_params:
            self.protein_params = {
                f"R{i + 1:04d}": {
                    "cx": 0.0,
                    "cy": 0.0,
                    "cz": 0.0,
                    "rx": 0.0,
                    "ry": 0.0,
                    "rz": 0.0,
                }
                for i in range(self.n_systems)
            }

    def to_json(self, path: str):
        """Save the settings to a JSON file, without the session-only values.

        Args:
            path: The file to write.
        """
        data = asdict(self)

        transient_keys = {
            "lipid_entries_absolute",
            "lipid_entries_relative",
            "membrane_string",
            "entries_current",
            "extern_config",
            "external_lipid_definitions",
        }
        for key in transient_keys:
            data.pop(key, None)

        with open(path, "w") as f:
            json.dump(data, f, indent=4)

    @classmethod
    def from_json(cls, path: str):
        """Read settings saved with `to_json`.

        Args:
            path: The JSON file.

        Returns:
            A new `Config`.
        """
        with open(path, "r") as f:
            data = json.load(f)
        return cls(**data)
