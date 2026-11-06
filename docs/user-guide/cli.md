# Command line (`omb cli`)

Every build writes a `config.json` with all its settings. `omb cli` builds systems from such a
file, without opening a browser:

```bash
omb cli output/OpenMembraneBuilder-20261106_101500/config.json
```

The result is a new output folder `output/OpenMembraneBuilder-<date>_<time>-headless/`, with
the same structure as a build from the app, its own copy of `config.json` and a zip archive in
`downloads/`.

## Why

- **Replicas** – change `n_systems` and build more systems with the same settings.
- **Variations** – change the box, the composition, the salt or the protein placement in the
  file, and build again.
- **Clusters and scripts** – no browser needed; run it where the simulations run.
- **Reproducibility** – keep the `config.json` with your data; it rebuilds the setup.

## Typical use

```bash
# 1. design the system once in the app
omb gui

# 2. rebuild / vary it from the file
cp output/OpenMembraneBuilder-20261106_101500/config.json popc_ubq.json
#    … edit popc_ubq.json (e.g. "n_systems": 10) …
omb cli popc_ubq.json
```

### Several variations in one go

```bash
for salt in 0.05 0.15 0.30; do
  sed "s/salt_molarity:[0-9.]*/salt_molarity:${salt}/" base.json > salt_${salt}.json
  omb cli salt_${salt}.json
done
```

### Elsewhere

```bash
OMB_WORKDIR=/scratch/$USER/builds omb cli base.json     # output/ and downloads/ go there
python -m openmembranebuilder cli base.json             # the same, if `omb` is not on the PATH
```

## What is rebuilt, and what is reused

The file has one section per step of the app (see [config.json](../reference/config.md)):

- **`martinize`** – how the proteins were prepared. It is **not** run again: the finished
  protein files are recorded in the `membrane` section (`protein_files`) and used directly.
- **`membrane`** – everything else; this is what is built.

!!! important "Run it from the same folder"
    The protein files in `protein_files` point into the earlier build (e.g.
    `output/OpenMembraneBuilder-…/user_inputs/1UBQ.pdb`), relative to the folder where that
    build ran. Start `omb cli` from the same folder, or edit the paths to where the files are
    now.

## What you can change

| To change … | edit |
|---|---|
| number of systems | `n_systems` |
| box | `box_type`, `box_x`, `box_y`, `box_z` |
| composition | `entries` (and `abs_lip_vals` for numbers of lipids) – see [config.json](../reference/config.md#the-composition) |
| several compositions at once | `entries` as `{"name": [...], ...}` – see [Membrane templates](templates.md#using-a-template) |
| salt and ions | `solvation`, e.g. `"solv:W pos:NA neg:CL salt_molarity:0.15"` |
| force field | `selected_ff` (`martini_v3` / `martini_v2.2`) – the proteins must match |
| protein files | `protein_files` |
| protein placement | `protein_params`, `z_method`, `distance_to_mem`, `randomize_*`, `copy_*` |

The placement used for every system is in `protein_params`, so a rebuild places the
proteins exactly as before. With random placements (`randomize_pos`, `randomize_rot`) new
random values are drawn.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | built |
| 2 | wrong command line (e.g. a missing file name) |
| 3 | the configuration file does not exist |
