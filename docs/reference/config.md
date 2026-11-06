# config.json

Every build writes `config.json` into its output folder (and into the zip). It records every
setting, so the same systems can be built again with [`omb cli`](../user-guide/cli.md), and it
documents how a system was made.

The file has one section per step of the app:

```text
{
    "martinize": { ... },     only with proteins: how they were prepared
    "membrane":  { ... }      the build: box, water, composition, proteins, placement
}
```

`omb cli` builds the `membrane` section. The `martinize` section is kept for the record; it is
not run again, because the finished protein files are listed in `membrane.protein_files`.

## A minimal file

A membrane without proteins needs only a few keys; everything else takes its default.

```json title="popc_chol.json"
{
    "membrane": {
        "selected_ff": "martini_v3",
        "box_type": "rectangular",
        "box_x": 12.0, "box_y": 12.0, "box_z": 12.0,
        "solvation": "solv:W pos:NA neg:CL salt_molarity:0.15",
        "n_systems": 3,
        "abs_lip_vals": false,
        "entries": [
            ["POPC", 0.7, 0.7, 0.64, 0.64],
            ["CHOL", 0.3, 0.3, 0.40, 0.40]
        ]
    }
}
```

```bash
omb cli popc_chol.json
```

## The `membrane` section

### General

| Key | Type | Default | Meaning |
|---|---|---|---|
| `selected_ff` | text | `"martini_v3"` | force field: `"martini_v3"` or `"martini_v2.2"` |
| `n_systems` | number | `1` | number of independent systems (`R0001`, `R0002`, …) |
| `output_name` | text | `""` | recorded name of the output; a rebuild always writes a new `…-headless` folder |

### Box

| Key | Type | Default | Meaning |
|---|---|---|---|
| `box_type` | text | `"rectangular"` | `"rectangular"`, `"hexagonal"` or `"skewed_hexagonal"` |
| `box_x` | nm | `10.0` | length x (all box types) |
| `box_y` | nm | `10.0` | length y (rectangular only) |
| `box_z` | nm | `10.0` | height z (rectangular and hexagonal) |

### Water and ions

| Key | Type | Default | Meaning |
|---|---|---|---|
| `solvation` | text | | COBY solvation text, see below |
| `salt_molarity` | M | `0.15` | recorded salt concentration (the value used is the one in `solvation`) |

```text
solv:W pos:NA neg:CL salt_molarity:0.15                                  Martini 3 / 2.2
solv:W90.0 solv:WF10.0 params:IMPORTED pos:NA neg:CL salt_molarity:0.15  Martini 2.2 with 10 % antifreeze water
```

Positive ions: `NA`, `CA`. Negative ions: `CL`, `BR`, `IOD`, `ACE`, `BF4`, `PF6`, `SCN`, `CLO4`,
`NO3`.

### The composition

| Key | Type | Default | Meaning |
|---|---|---|---|
| `entries` | list or object | `[]` | the lipids, see below |
| `abs_lip_vals` | true/false | `false` | `false`: leaflet values are fractions (each leaflet adds up to 1); `true`: numbers of lipids |
| `template_active` | true/false | `false` | the composition came from a template |
| `template_path` | text | `""` | that template (copied to `user_inputs/`) |

Each lipid is one row:

```text
[ lipid,  upper,  lower,  area per lipid upper,  area per lipid lower ]
[ "POPC", 0.7,    0.7,    0.64,                  0.64 ]
```

Areas per lipid are in nm². `entries` can also hold **several compositions**, each built with
its own systems (`systems/<name>_R0001/`, …):

```json
"entries": {
    "PlasmaMembrane":  [["POPC", 0.7, 0.3, 0.61, 0.61], ["DOPC", 0.3, 0.7, 0.6, 0.6]],
    "CholesterolRich": [["POPC", 0.5, 0.4, 0.61, 0.61], ["CHOL", 0.5, 0.6, 0.4, 0.4]]
}
```

### Proteins

| Key | Type | Meaning |
|---|---|---|
| `protein_files` | object | the proteins: `{"1UBQ.pdb": {"pdb_path": …, "itp_path": …, "force_field": "martini_v3"}}`. The key is the protein's file name; its name (without `.pdb`) is the molecule name in the topology. |

The paths are relative to the folder where the build ran (or absolute).

### Protein placement

| Key | Type | Default | Meaning |
|---|---|---|---|
| `z_method` | text | `"Absolute z position"` | or `"Height above Membrane"` |
| `distance_to_mem` | nm | `2.0` | gap between membrane surface and lowest protein bead (*Height above Membrane*) |
| `randomize_pos` | true/false | `false` | random x/y position |
| `randomize_pos_every` | true/false | `false` | a new random position for every system |
| `randomize_rot` | true/false | `false` | random orientation |
| `randomize_rot_every` | true/false | `false` | a new random orientation for every system |
| `use_protein_distance` | true/false | `false` | two proteins: place them by their gap |
| `protein_distance` | nm | `1.0` | that gap |
| `more_copies` | true/false | `false` | several copies of one protein |
| `copy_number` | number | `2` | copies per system |
| `copy_mode` | text | `"together"` | `"together"` (grid, shared orientation) or `"independent"` (own pose per copy) |
| `copy_params` | list | `[]` | independent copies: `[{"cx", "cy", "cz", "rx", "ry", "rz"}, …]` |
| `randomize_rot_all_copies` | true/false | `false` | copies "together": a random orientation for each copy |
| `protein_params` | object | | the pose of every protein in every system, see below |
| `insertion_params` | object | | the settings of the protein insertion page, per protein (used to restore the page) |

`protein_params` holds the position (nm from the box centre) and rotation (degrees) used for
each protein in each system:

```json
"protein_params": {
    "R0001": {"1UBQ.pdb": {"cx": 0.0, "cy": 0.0, "cz": 4.12, "rx": 0.0, "ry": 0.0, "rz": 0.0}},
    "R0002": {"1UBQ.pdb": {"cx": 0.0, "cy": 0.0, "cz": 4.12, "rx": 0.0, "ry": 0.0, "rz": 0.0}}
}
```

Independent copies add a `"copies": [ … ]` list with the pose of every copy.
Fixed (not random) positions and rotations are taken from `R0001` for all systems.

### Recorded for the app

| Key | Meaning |
|---|---|
| `molecule_builder_index` | custom lipids made by a molecule builder (empty unless used) |
| `extern_config` | settings handed in by another tool (empty unless used) |

## The `martinize` section

One entry per protein in `proteins`, with how it was prepared. `omb cli` uses it to prepare
the proteins again with `--martinize` or `--pdb`, or when `membrane.protein_files` is empty
(see [Command line](../user-guide/cli.md#what-is-rebuilt-and-what-is-reused)):

| Key | Meaning |
|---|---|
| `name` | the protein's name (`1UBQ`) |
| `protein_input_mode` | `"upload_cg"` or `"martinize"` |
| `cg_pdb_path`, `cg_itp_path`, `cg_martini_version` | uploaded coarse-grained files and their Martini version |
| `atomistic_source`, `atomistic_pdb_path`, `fetch_pdb_id` | where the atomistic structure came from |
| `clean_structure`, `selected_chains`, `chain_ranges` | clean-up, chains and residue ranges kept |
| `do_mutation`, `mutation_chain`, `mutation_residue_idx`, `mutation_target` | the mutation, if any |
| `forcefield`, `secondary_structure_mode`, `custom_ss_string` | martinize2 force field and secondary structure |
| `use_elastic_network`, `elastic_lower`, `elastic_upper`, `elastic_force` | elastic network |
| `position_restraints`, `cysteine_bridges`, `extra_flags` | further martinize2 options |
| `output_cg_pdb_path`, `output_itp_path` | the finished `<name>.pdb` / `<name>.itp` |

## A complete example

??? example "config.json of a build with two independent copies of 1UBQ"

    ```json
    --8<-- "docs/reference/example_config.json"
    ```

## Notes

- Values that only matter during a session (temporary folders, intermediate text) are not
  written.
- The `membrane` section accepts only the keys listed on this page; an unknown key stops
  the build with an error naming it. (Old keys in the `martinize` section are ignored.)
- Keep the file with your data: together with the program version (`omb --version`) it
  describes exactly how the systems were made.
