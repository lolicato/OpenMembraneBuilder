# Installation

OpenMembraneBuilder is a Python package. Installing it gives you the `omb` command:

```bash
omb gui                 # the graphical app
omb cli config.json     # a build from a saved configuration
```

---

## 📦 From the compressed file

```bash
pip install openmembranebuilder-1.2.0.dev0.tar.gz
# or the wheel:
pip install openmembranebuilder-1.2.0.dev0-py3-none-any.whl
```

All Python dependencies (COBY, MDAnalysis, Streamlit, vermouth/martinize2, …) are installed
with it. Python 3.10 or newer is required.

---

## 📥 From the repository

### Stable version

```bash
git clone https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
pip install .
```

### Development version (specific branch)

```bash
git clone -b development https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
pip install -e .        # editable: code changes take effect without reinstalling
```

---

## 🐍 Conda Environment

A Conda environment is still a good home for it:

```bash
conda env create -f environment.yml     # dependencies + OpenMembraneBuilder itself
conda activate OpenMembraneBuilder      # OpenMembraneBuilder-dev for the development version
```

---

## 🚀 Run OpenMembraneBuilder

### GUI mode

```bash
omb gui
omb gui --port 8502 --server.headless true   # extra options go to Streamlit
```

### CLI mode

```bash
omb cli outputs/OpenMembraneBuilder-xxxx/config.json
```

Builds are written to the folder you start `omb` from. Set `OMB_WORKDIR` to write them
somewhere else:

```bash
OMB_WORKDIR=/scratch/me/membranes omb cli config.json
```

The old way, `streamlit run app.py` and `python app.py --no-gui config.json` from the
repository folder, still works.

---

## 📦 Notes

- `coby` is used for system generation
- `MDAnalysis` handles structure manipulation
- `streamlit` provides the GUI
- `streamlit-molstar` enables 3D visualization
- `vermouth` provides `martinize2` for coarse-graining proteins
- PyMOL (only for mutations) is installed separately:
  `conda install -c conda-forge pymol-open-source`

---

## 🔄 Updating

```bash
pip install --upgrade openmembranebuilder-X.Y.Z.tar.gz
```

---

## 🔧 Optional (Documentation)

If you want to build the documentation locally:

```bash
pip install ".[docs]"
mkdocs serve
```
