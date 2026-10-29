from dataclasses import dataclass, asdict, field
import json


@dataclass
class Config:
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
        with open(path, "r") as f:
            data = json.load(f)
        return cls(**data)
