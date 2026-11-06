# Installation

OpenMembraneBuilder is a Python package. Installing it gives you one command, `omb`:

```bash
omb gui                 # the graphical app, in your browser
omb gui --no-run        # the same, without the button that starts the simulations
omb cli config.json     # a build from a saved configuration, without the browser
omb --version           # the installed version
```

## Requirements

| | |
|---|---|
| Python | 3.10 or newer |
| Operating system | Linux, macOS or Windows (WSL recommended for running GROMACS) |
| Python packages | installed automatically: COBY, MDAnalysis, mdtraj, vermouth (martinize2), Streamlit, streamlit-molstar, NumPy, SciPy, pandas |
| Optional | **PyMOL** – only to mutate residues; **jq** – used by `run.sh` to create the index file; **GROMACS** – to run the simulations |

## Install

=== "From the compressed file"

    ```bash
    pip install openmembranebuilder-2.0.0.tar.gz
    # or the wheel – the same package, already built:
    pip install openmembranebuilder-2.0.0-py3-none-any.whl
    ```

=== "From the repository"

    ```bash
    git clone -b development https://github.com/lolicato/OpenMembraneBuilder.git
    cd OpenMembraneBuilder
    pip install .          # or: pip install -e .   (editable: code changes apply at once)
    ```

=== "With conda"

    ```bash
    git clone -b development https://github.com/lolicato/OpenMembraneBuilder.git
    cd OpenMembraneBuilder
    conda env create -f environment.yml     # dependencies + OpenMembraneBuilder itself
    conda activate OpenMembraneBuilder
    ```

!!! tip "A separate environment"
    Installing into its own environment keeps OpenMembraneBuilder's dependencies away from
    other projects:

    ```bash
    python -m venv omb-env && source omb-env/bin/activate     # or: conda create -n omb python=3.12
    pip install openmembranebuilder-2.0.0.tar.gz
    ```

### Optional tools

```bash
conda install -c conda-forge pymol-open-source   # residue mutations
conda install -c conda-forge jq                   # or: brew install jq / sudo apt install jq
```

GROMACS is only needed where the simulations run (for example on a cluster); building does not
need it.

### Check

```bash
omb --version            # OpenMembraneBuilder 2.0.0
martinize2 -h | head -3  # comes with the vermouth package
```

## Where your files go

`omb` writes everything into the folder you start it from:

| Folder | Content |
|---|---|
| `output/` | one sub-folder per build, with all systems |
| `downloads/` | a zip archive of every build |
| `temp-uploads_…/` | uploads during a session; removed after the build |

To write somewhere else, set `OMB_WORKDIR`:

```bash
OMB_WORKDIR=/scratch/me/membranes omb gui
```

The force fields, simulation settings, run scripts and membrane templates come with the
package; nothing has to be downloaded.

## Updating

A file with a higher version number replaces the installed one with a plain
`pip install`. To reinstall a rebuilt file with the same number:

```bash
pip install --force-reinstall --no-deps openmembranebuilder-2.0.0.tar.gz
```

## Making the compressed file yourself

From the repository folder:

```bash
pip install build
python -m build      # writes dist/openmembranebuilder-2.0.0.tar.gz and the .whl
```

The `.tar.gz` contains the source; the `.whl` (a "wheel") is the same package already built,
so it installs a little faster. Both work on every operating system.

## Building this documentation

```bash
pip install ".[docs]"
mkdocs serve        # live preview on http://127.0.0.1:8000
mkdocs build        # static site in site/
```

!!! note "The previous way still works"
    From the repository folder, `streamlit run app.py` and `python app.py --no-gui config.json`
    behave as before.
