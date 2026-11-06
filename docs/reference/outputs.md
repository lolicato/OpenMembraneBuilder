# Output and GROMACS

## The output folder

Each build creates one folder in `output/` and a zip of it in `downloads/`, in the folder
where `omb` was started (or in `$OMB_WORKDIR`).

```text
output/OpenMembraneBuilder-20261106_101500-ubq/
├── config.json              every setting of this build
├── user_inputs/             your input files, as provided/prepared
│   ├── 1UBQ.pdb
│   ├── 1UBQ.itp
│   └── my_template.ob       (when a template was uploaded)
├── toppar/                  everything the topology includes
│   ├── martini_v3/          the force field (.itp files)
│   └── 1UBQ.itp             the protein topology ([ moleculetype ] 1UBQ 1)
├── mdp/                     simulation settings
│   ├── em.mdp  eq1.mdp  eq2.mdp  eq3.mdp  md.mdp
├── scripts/
│   └── index.sh             makes the index groups (used by run.sh)
└── systems/
    ├── run.sh               runs all systems
    ├── run.log, run.pid     only after "Run simulations" in the app
    ├── R0001/
    │   ├── system.gro       the built system
    │   ├── system.pdb       the same, for viewers
    │   ├── topol.top        its topology
    │   └── 1UBQ.pdb         the protein file COBY read
    ├── R0002/ …
downloads/OpenMembraneBuilder-20261106_101500-ubq.zip
```

| Folder name | When |
|---|---|
| `OpenMembraneBuilder-<date>_<time>` | built in the app |
| `OpenMembraneBuilder-<date>_<time>-<name>` | built in the app with an output folder name |
| `OpenMembraneBuilder-<date>_<time>-headless` | built with `omb cli` |
| `OpenMembraneBuilder-<date>_<time>-headless-1UBQ` | built with `omb cli … --pdb 1UBQ` |

A build from `omb cli` that prepared its proteins also contains `martinize/p0/`, … with the
work of the protein preparation (structures, `martinize.log`, the martinize2 command).

System folders are `R0001`, `R0002`, …; with several compositions (see
[Membrane templates](../user-guide/templates.md#using-a-template)) they are
`<composition>_R0001`, ….

## The topology

```text title="systems/R0001/topol.top"
#include "../../toppar/martini_v3/martini_v3.0.0.itp"
#include "../../toppar/martini_v3/martini_v3.0.0_ffbonded_v2.itp"
#include "../../toppar/martini_v3/martini_v3.0.0_solvents_v1.itp"
#include "../../toppar/martini_v3/martini_v3.0.0_ions_v1.itp"
#include "../../toppar/martini_v3/martini_v3.0.0_phospholipids_PC_v2.itp"
...
#include "../../toppar/1UBQ.itp"

[ system ]
1UBQ_POPC_0.7_0.7_CHOL_0.3_0.3

[ molecules ]
; name number
1UBQ        2
POPC        449
CHOL        192
W           13171
NA          143
CL          143
```

- **Includes** point into `toppar/` of the same build, so the folder can be moved or copied to
  a cluster as a whole.
- **Title**: the protein names and the composition (lipid, upper, lower for each lipid).
- **Proteins** appear under their own name; copies of one protein are included once and
  counted in one line.
- **Order** of the molecules is the order in `system.gro`.

## Index groups

`run.sh` creates `index.ndx` after the minimisation (with `scripts/index.sh` and
`gmx make_ndx`). The groups are used by the thermostat (`tc-grps = SOLU SOLV`).

| Group | Content |
|---|---|
| `MEMB` | all lipids |
| `SOLV` | water and ions |
| `PROT` | the proteins (everything that is neither membrane nor solvent) – only with proteins |
| `SOLU` | proteins + membrane (or the membrane alone, without proteins) |

## The simulation protocol

The `mdp/` folder holds the protocol for the chosen force field; systems with proteins get
the `protein` set, membranes the `membrane` set (same steps).

| Step | File | Method | Time step | Steps | Length | Pressure |
|---|---|---|---|---|---|---|
| Minimisation | `em.mdp` | steepest descent, flexible water | – | 10 000 | – | – |
| Equilibration 1 | `eq1.mdp` | MD | 5 fs | 50 000 | 0.25 ns | none (NVT) |
| Equilibration 2 | `eq2.mdp` | MD | 10 fs | 50 000 | 0.5 ns | C-rescale, semi-isotropic, 1 bar |
| Equilibration 3 | `eq3.mdp` | MD | 20 fs | 100 000 | 2 ns | C-rescale, semi-isotropic, 1 bar |
| Production | `md.mdp` | MD | 20 fs | 250 000 000 | 5 µs | C-rescale, semi-isotropic, 1 bar |

All steps: 310 K (v-rescale on `SOLU` and `SOLV`), reaction-field electrostatics with
ε<sub>r</sub> = 15, 1.1 nm cut-offs, Verlet lists – the standard Martini settings. Adjust the
files (e.g. `nsteps` of `md.mdp`, the temperature) before running.

## Running

From the app: press **▶ Run simulations** on the results page
([details](../user-guide/build-results.md#run-the-simulations)). By hand:

```bash
cd output/OpenMembraneBuilder-20261106_101500-ubq/systems
./run.sh
```

`run.sh` goes through every system folder and, in each one, runs

1. `gmx grompp` + `gmx mdrun` for **em** (from `system.gro`),
2. creates **index.ndx** (needs `jq`),
3. **eq1**, **eq2**, **eq3** (each starting from the previous step, with `-r` as reference for
   restraints) and **md**.

Each step writes `<step>.tpr`, `.log`, `.edr`, `.gro` … into the system folder.

!!! tip "Before you run"
    - `gmx` must be on the PATH (`source /path/to/gromacs/bin/GMXRC`).
    - `jq` must be installed (`conda install -c conda-forge jq`).
    - Start `run.sh` from the `systems/` folder: it runs every sub-folder of the current folder.
    - On a cluster, wrap it in your batch script, or run the `gmx` commands per system.

!!! note "“Energy minimization reached the maximum number of steps”"
    With proteins, the minimisation often stops before the forces fall below the (strict)
    default target of 10 kJ mol⁻¹ nm⁻¹. This is normal for Martini; GROMACS still writes the
    lowest-energy structure and `run.sh` continues with it. To silence the message, add
    `emtol = 1000` to `em.mdp`.

## The zip archive

`downloads/<folder>.zip` contains the whole output folder, including `config.json`. Unpack it
anywhere – all paths inside are relative.
