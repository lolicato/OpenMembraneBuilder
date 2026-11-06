# Force fields and settings

Everything a build needs is shipped inside the package, in `openmembranebuilder/resources/`:

```text
resources/
├── toppar/                     force fields
│   ├── martini_v3.top          which .itp files make up Martini 3
│   ├── martini_v3/             Martini 3.0.0 (lipids, sterols, solvents, ions, bonded)
│   ├── martini_v2.2.top
│   └── martini_v2.2/           Martini 2.2 (lipids, cholesterol, solvents, ions)
├── mdp/<force field>/
│   ├── membrane/               em, eq1–eq3, md for membranes
│   └── protein/                the same for systems with proteins
├── scripts/                    run.sh, index.sh
├── membrane_templates/         predefined compositions (.csv / .ob)
├── example_files/              example_template.ob (offered for download in the app)
└── file_imports/martini_v2.2/  CHOL.pdb, CHOL2.pdb, WF.gro (Martini 2.2 structures)
```

## The force fields

| | `martini_v3` | `martini_v2.2` |
|---|---|---|
| Lipids | PC, PE, PG, PS, PA, PI and phosphoinositides, SM, cardiolipin, sterols | PC, PE, PG, PS, PA, SM, cholesterol |
| Water | `W` | `W`, optionally part antifreeze `WF` |
| Proteins (martinize2) | `martini3001`, `martini3001-idp` | `martini22` |

The lipid list in the app comes from these files: every `.itp` included by the `.top` file
whose name points to lipids (phospholipids, sterols, …) is read, and each molecule it defines
is offered. **Show Lipid Mapping** on the membrane page describes each lipid (type, head
group, linker, tails, bead types), taken from the comments in the files.

### Martini 2.2 specifics

- **Cholesterol**: `CHOL` (the default model) or `CHOL2` (an improved model); its structure
  is provided by the program.
- **Antifreeze water**: *Percentage of unfreezable water* replaces part of `W` by `WF`.

## Simulation settings

See [Output and GROMACS](outputs.md#the-simulation-protocol) for the protocol. To change it
for all future builds, edit the files in `resources/mdp/` of your installation (or of your
source checkout before `pip install`).

## Adding lipids or a force field

A force field is a `<name>.top` file in `resources/toppar/` that includes `.itp` files from
`resources/toppar/<name>/`, plus a matching `resources/mdp/<name>/membrane|protein/` folder.
New lipids are added by adding their `.itp` (with a lipid-type word such as `phospholipids`
or `sterols` in the file name) and including it in the `.top` file.

!!! note
    COBY must know how to build a lipid. Lipids from the Martini 3 / 2.2 libraries are known
    to COBY; for others see the COBY documentation.

## Checking what is installed

```bash
python -c "import openmembranebuilder, os; print(os.path.dirname(openmembranebuilder.__file__))"
```

prints the package folder; `resources/` is inside it.
