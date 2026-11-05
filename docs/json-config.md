# JSON Configuration

Each build writes a `config.json` file inside the output folder.

This file stores the full build configuration and allows full reproducibility of the system.

---

## 📄 Example Configuration

```json
{
    "output_name": "outputs/OpenMembraneBuilder-2026-05-15_00-47-52",
    "selected_ff": "martini_v3",
    "base_folder": "",
    "selected_module": "membrane_with_cg_protein",
    "box_x": 10.0,
    "box_y": 10.0,
    "box_z": 20.0,
    "box_type": "rectangular",
    "solvation": "solv:W pos:NA neg:CL salt_molarity:0.18",
    "salt_molarity": 0.18,
    "n_systems": 4,
    "pdb_path": "outputs/OpenMembraneBuilder-2026-05-15_00-47-52/simulations/R0004/protein.pdb",
    "itp_path": "outputs/OpenMembraneBuilder-2026-05-15_00-47-52/toppar/protein.itp",
    "abs_lip_vals": false,
    "entries": [
        [
            "POPC",
            0.3,
            0.25,
            0.6,
            0.6
        ],
        [
            "CHOL",
            0.3,
            0.5,
            0.6,
            0.6
        ],
        [
            "PSM",
            0.4,
            0.25,
            0.6,
            0.6
        ]
    ],
    "z_method": "Height above Membrane",
    "distance_to_mem": 2.0,
    "randomize_pos": true,
    "randomize_pos_every": true,
    "randomize_rot": true,
    "randomize_rot_every": true,
    "protein_params": {
        "R0001": {
            "cx": -0.2634,
            "cy": -0.3105,
            "cz": 5.5903,
            "rx": 140.6176,
            "ry": 154.579,
            "rz": -117.2219
        },
        "R0002": {
            "cx": 1.7503,
            "cy": 0.6474,
            "cz": 5.9287,
            "rx": 72.6335,
            "ry": -121.346,
            "rz": -131.9327
        },
        "R0003": {
            "cx": -1.6669,
            "cy": 1.3716,
            "cz": 5.8229,
            "rx": 17.8925,
            "ry": -115.151,
            "rz": 177.4764
        },
        "R0004": {
            "cx": -0.9433,
            "cy": -1.1002,
            "cz": 5.914,
            "rx": -169.1105,
            "ry": 45.0814,
            "rz": -55.516
        }
    }
}
```

---

## 🧠 Explanation of Key Fields

### 🧱 System

- `selected_module`  
  - `"membrane"`  
  - `"membrane_with_cg_protein"`

- `selected_ff`  
  - Force field (currently `martini_v3`)

---

### 📦 Box

- `box_x`, `box_y`, `box_z` → size in nm  
- `box_type` → currently `"rectangular"`

---

### 🧂 Solvation

- `solvation` → string passed to COBY  
- `salt_molarity` → numeric version (used internally)

---

### 🔁 Replicas

- `n_systems` → number of independent systems generated

---

### 🧬 Protein

- `pdb_path` → protein structure  
- `itp_path` → topology  

These are automatically copied into:

```
user_inputs/
```

---

### 📍 Placement

- `z_method`
  - `"Absolute z position"`
  - `"Height above Membrane"`

- `distance_to_mem`
  - Used only for height-based placement

---

### 🎲 Randomization

- `randomize_pos`
- `randomize_pos_every`
- `randomize_rot`
- `randomize_rot_every`

Controls variability across replicas.

---

### 🧬 Membrane Composition

Each lipid entry follows:

```
[lipid_name, upper, lower, apl_upper, apl_lower]
```

Example:

```json
["POPC", 0.5, 1.0, 0.6, 0.6]
```

#### Relative mode (`abs_lip_vals = false`)
- Values must sum to **1.0 per leaflet**

#### Absolute mode (`abs_lip_vals = true`)
- Values are **lipid counts (integers)**

---

## ⚠️ Notes

- `membrane_string` is generated automatically → do not edit  
- `entries` is an internal field  
- `output_name` is set automatically  

---

## 🚀 Reuse

```bash
omb cli outputs/OpenMembraneBuilder-xxxx/config.json
```
