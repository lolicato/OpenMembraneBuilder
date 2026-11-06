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

- **`martinize`** – how the proteins are prepared: where the structure comes from, clean-up,
  chains, mutation and the martinize2 options.
- **`membrane`** – the system: box, water, lipids, the finished protein files
  (`protein_files`) and the placement of the proteins.

By default the finished protein files recorded in `protein_files` are **reused** – martinize2
is not run again, and the systems are built exactly from the proteins you prepared before.

The proteins are **prepared again** (download or file → clean-up → chains → mutation →
martinize2, exactly as in the app) when

| | |
|---|---|
| `--pdb ID_OR_FILE` | the protein is replaced by another structure (see below) |
| `--martinize` | you ask for it, e.g. after changing martinize2 options in the file |
| no `protein_files` | the file has a `martinize` section but no finished proteins (e.g. a template you wrote by hand) |

The work of the preparation is kept in the new output folder, under `martinize/p0/`,
`martinize/p1/`, … (structures, `martinize.log`, the exact martinize2 command), and the
`config.json` written there records the prepared files – so that output rebuilds again
without running martinize2.

!!! note "Proteins uploaded as atomistic files in the app"
    The app keeps uploaded atomistic structures only during the session. To prepare such a
    protein again from the command line, pass the file: `--pdb my_protein.pdb`. (Proteins
    fetched by PDB code are simply downloaded again.)

!!! important "Where the protein files are looked for"
    `protein_files` points into the earlier build (`…/user_inputs/1UBQ.pdb`, `…/toppar/1UBQ.itp`).
    If you moved or deleted that build, use `--martinize` (or `--pdb`) to prepare the protein
    again, or edit the paths.

## Many proteins with one template

`--pdb` turns any configuration into a template: the protein is replaced by another structure,
which is prepared with the martinize2 settings of the file, and inserted with its placement
(height above the membrane, rotation, copies …).

```bash
# 1. set everything up once in the app (e.g. with 1UBQ) and keep its config.json
cp output/OpenMembraneBuilder-20261106_101500/config.json template.json

# 2. the same system for other proteins
for PDB in 1BFF 1UBQ 2RH1 4HHB; do
    omb cli template.json --pdb ${PDB}
done
```

Each run creates its own folder, named after the protein:
`output/OpenMembraneBuilder-<date>_<time>-headless-1BFF/`, … .

| `--pdb` value | Source |
|---|---|
| a PDB code, e.g. `1UBQ` | downloaded from the RCSB PDB |
| a file, e.g. `structures/fgf2.pdb` | your atomistic structure |

- The protein's **name** becomes the PDB code (upper case) or the file name without `.pdb`
  (`fgf2`), so the files and the topology are called `1UBQ.pdb`, `1UBQ.itp`, `1UBQ`.
- **Kept from the template**: clean-up, martinize2 force field, secondary structure, elastic
  network, restraints, disulfide bridges, extra options – and the protein's placement.
- **Reset**: chain and residue selections and the mutation, which belong to one particular
  structure – all chains of the new structure are used.
- **Several proteins**: give `--pdb` once per protein, in order
  (`--pdb 1BFF --pdb 1UBQ`); extra ones get the settings of the last protein.

Local files work the same way:

```bash
for f in structures/*.pdb; do
    omb cli template.json --pdb "$f"
done
```

!!! tip "When a protein fails"
    `omb cli` stops with `Error: <name>: …` and exit code 1 (e.g. a code that does not exist,
    no internet access, or martinize2 failing on a structure with gaps), and the loop goes on
    with the next protein. Collect the failures with
    `omb cli template.json --pdb $PDB || echo $PDB >> failed.txt`.

## What you can change

| To change … | edit |
|---|---|
| number of systems | `n_systems` |
| box | `box_type`, `box_x`, `box_y`, `box_z` |
| composition | `entries` (and `abs_lip_vals` for numbers of lipids) – see [config.json](../reference/config.md#the-composition) |
| several compositions at once | `entries` as `{"name": [...], ...}` – see [Membrane templates](templates.md#using-a-template) |
| salt and ions | `solvation`, e.g. `"solv:W pos:NA neg:CL salt_molarity:0.15"` |
| force field | `selected_ff` (`martini_v3` / `martini_v2.2`) – the proteins must match |
| protein files | `protein_files` – or another structure with `--pdb` |
| martinize2 options | `martinize` → `proteins` (then use `--martinize`) |
| protein placement | `protein_params`, `z_method`, `distance_to_mem`, `randomize_*`, `copy_*` |

The placement used for every system is in `protein_params`, so a rebuild places the
proteins exactly as before. With random placements (`randomize_pos`, `randomize_rot`) new
random values are drawn.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | built |
| 1 | a protein could not be prepared (the message says why) |
| 2 | wrong command line (e.g. a missing file name) |
| 3 | the configuration file does not exist |
