# Troubleshooting

## Installing and starting

??? question "`omb: command not found`"
    The environment in which you installed OpenMembraneBuilder is not active, or its `bin`
    folder is not on the PATH. Activate it (`conda activate …` / `source …/bin/activate`), or
    run `python -m openmembranebuilder gui`.

??? question "The browser does not open / port 8501 is in use"
    Open the address printed in the terminal by hand. If the port is taken, use another one:
    `omb gui --port 8502`. On a remote machine use `omb gui --server.headless true` and an SSH
    tunnel (see [The workflow](../user-guide/workflow.md#starting-the-app)).

??? question "Where are my results?"
    In `output/` and `downloads/` of the folder where you started `omb`, or in `$OMB_WORKDIR`
    if it is set.

## Proteins

??? question "Name: use only letters, digits, _ and -"
    Names become file names and topology names; spaces and dots are not allowed. `GSDMD_NT`
    or `GSDMD-NT` work, `GSDMD NT` and `gsdmd.1` do not.

??? question "Each protein needs its own name"
    Two proteins have the same name (upper/lower case does not count). Give the second one a
    different name, e.g. `FGF2_A` and `FGF2_B`.

??? question "martinize2 failed"
    The results page shows the command and martinize2's message. Common causes:

    - missing atoms or non-standard residues – clean the structure, cut away the problematic
      residues with the residue range, or repair the structure first;
    - a custom secondary structure with the wrong length – it needs one letter per residue;
    - an unknown option in *Extra martinize2 options*.

??? question "Mismatch found between number of beads in structure and number of beads in topology"
    martinize2 made **more than one molecule** out of the protein – because it has several
    chains, or gaps (missing residues) that cut a chain into pieces – but only the topology of
    the first piece is handed on, so the structure has more beads than the topology.

    What helps:

    - keep one chain without gaps: choose a single chain and a residue range between the gaps;
    - fill the gaps first (e.g. with PDBFixer or MODELLER) and use the repaired structure;
    - ask martinize2 to make one molecule of the chains with `-merge A` (or `-merge A,B`) in
      *Extra martinize2 options*;
    - or run martinize2 yourself and upload the result as an already coarse-grained protein.

??? question "No such file or directory: 'martinize2'"
    martinize2 is installed with OpenMembraneBuilder, in the same environment. Start
    `omb gui` from that environment (activated), so its `bin` folder – which holds
    `martinize2` – is on the PATH.

??? question "PyMOL mutagenesis failed"
    PyMOL must be installed as the `pymol` command:
    `conda install -c conda-forge pymol-open-source`.

??? question "Could not fetch PDB ID"
    The code must have four characters and exist at RCSB, and the computer needs internet
    access. Otherwise download the file yourself and use *Upload PDB*.

## Membrane and placement

??? question "Ratios must sum to 1.0"
    With fractions, each leaflet column must add up to 1 (e.g. 0.7 + 0.3). For exact numbers
    of lipids tick *Give absolute lipid values*.

??? question "Lipid '…' is not available"
    The template names a lipid that the chosen force field does not have. Check the name in
    **Show Lipid Mapping**, or switch the force field.

??? question "A predefined template switched to absolute values"
    A template whose first line contains only whole numbers (also `1.0`) is read as numbers of
    lipids. Untick *Give absolute lipid values* after loading it.

??? question "… stick(s) out of the box"
    Part of a protein lies outside the box. Enlarge the box (on the same page), or move or
    rotate the protein. For *Height above Membrane*, check that Z leaves room above the
    membrane for the protein.

??? question "Overlapping: …"
    Two proteins (or copies) are closer than 0.3 nm somewhere. Move or rotate them apart, or
    enlarge the box.

??? question "Martini 2 / Martini 3 warning"
    The proteins were prepared for another Martini generation than the membrane force field.
    Prepare them again with the matching martinize2 force field (`martini3001` for
    `martini_v3`, `martini22` for `martini_v2.2`).

## Building and simulating

??? question "`omb cli`: unexpected keyword argument '…'"
    The `membrane` section contains a key this version does not know (often from an older
    file). Remove it; the valid keys are listed in [config.json](config.md).

??? question "`omb cli`: protein file not found"
    The paths in `protein_files` are relative to the folder where the original build ran.
    Start `omb cli` from that folder, or change the paths.

??? question "`jq: error … Cannot iterate over null` in run.sh"
    Fixed in 2.0.0: `index.sh` reads the current `config.json` format. In an older output,
    replace `scripts/index.sh` with the one from the package
    (`openmembranebuilder/resources/scripts/index.sh`).

??? question "`jq: command not found`"
    Install it: `conda install -c conda-forge jq`, `brew install jq` or `sudo apt install jq`.

??? question "Energy minimisation did not converge to Fmax < 10"
    Normal with proteins in Martini; GROMACS keeps the lowest-energy structure and `run.sh`
    continues. Add `emtol = 1000` to `mdp/em.mdp` to avoid the message
    (see [Output and GROMACS](outputs.md#running)).

??? question "The protein looks cut in two in VMD/PyMOL"
    It crosses the periodic box edge – nothing is wrong with the system. The results page
    of the app shows it in one piece; in GROMACS use `gmx trjconv -pbc mol`.

??? question "The Run simulations button is greyed out"
    `gmx` or `jq` was not found by the app. Start `omb gui` from a terminal in which both
    work (`gmx --version`, `jq --version`), e.g. after `source …/GMXRC`.

??? question "There is no Run simulations button"
    The app was started with `omb gui --no-run` (or `OMB_NO_RUN=1` is set) – on purpose, for
    shared servers. Start it without, or run `systems/run.sh` yourself.

## Reporting a problem

Please include the `config.json`, the output of `omb --version` and the error message (from
the page or the terminal) when you [open an issue](https://github.com/lolicato/OpenMembraneBuilder/issues).
