# Changelog

## 2.0.0

### Installation and commands

- Installable with `pip` – from the repository (`pip install .`) or from a single file
  (`openmembranebuilder-2.0.0.tar.gz` / `.whl`); all dependencies, including martinize2
  (vermouth), are installed with it.
- One command for everything: `omb gui`, `omb cli config.json`, `omb --version`
  (also `python -m openmembranebuilder …`).
- Builds are written to the folder `omb` is started from, or to `$OMB_WORKDIR`.
- The code lives in the `openmembranebuilder` package (formerly `src/`); force fields, mdp
  files, scripts and templates are installed with it. `streamlit run app.py` and
  `python app.py --no-gui config.json` still work.

### Run the simulations from the app

- **▶ Run simulations** on the results page starts `systems/run.sh` in the background after a
  build, with its log, a refresh and a stop button.
- `omb gui --no-run` (or `OMB_NO_RUN=1`) hides it, for apps offered on web servers.

### Protein insertion

- **System view**: box, membrane and proteins in x–z, y–z and x–y, drawn from the real
  protein beads and rotated exactly as COBY places them; warnings when a protein sticks out
  of the box or proteins overlap.
- The **box** can be changed on the protein insertion page.
- **Copies** of a protein can move together (grid) or **independently**, each with its own
  position and rotation.

### Proteins and topology

- Every protein gets a **name** (e.g. `1UBQ`): its files are `<name>.pdb` / `<name>.itp`, and
  the topology lists it as `<name>` (`[ moleculetype ]` and `[ molecules ]`).
- Copies of one protein are included **once** and counted in **one line** of `[ molecules ]`.

### Fixes

- `config.json` is included in the download zip; `omb cli` also writes it to the output.
- `index.sh` reads the current `config.json` format, so `run.sh` creates the index groups.
- Protein rotations follow COBY's convention (x, then y, then z about the turned axes).

### Documentation

- Every module, class and function explained in plain language; this documentation site
  with user guide, references and code reference.

## 1.0.3

- Command-line builds from a saved `config.json`.
- `config.json` written with every build and included in the zip archive.
- Temporary upload folders removed after each build.
- Cardiolipin lipids (TMCL, TOCL).
- Self-contained output folders with copies of the input files (`user_inputs/`).
