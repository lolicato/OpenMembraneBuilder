# 🚀 OpenMembraneBuilder-v2.0.0

![Version](https://img.shields.io/badge/version-v2.0.0-green)
![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)
![Status](https://img.shields.io/badge/status-active-success)
![License](https://img.shields.io/badge/license-AGPL--3.0-lightgrey)

👉 Working Progress!


🌐 **📖 Documentation:** https://lolicato.github.io/OpenMembraneBuilder  
👉 *Full user guide, and API reference available online*

---

**OpenMembraneBuilder** is a modern **GUI + CLI toolkit for building MARTINI coarse-grained membrane and membrane–protein systems**, powered by COBY.

It combines:
- 🧬 **MDAnalysis** for structure handling  
- ⚙️ **COBY** for system generation  
- 🧩 **martinize2** (vermouth) for coarse-graining proteins  
- 🌐 **Streamlit** for an interactive GUI  

with a strong focus on **reproducibility**, **usability**, and **local execution**.

---

## 🖼️ Interface

![Application GUI](./pictures/GUI.png)

---

## ✨ Features

- 🧪 Build membrane and membrane–protein systems  
- 🎛️ Fully interactive Streamlit GUI  
- 💻 One command for everything: `omb gui` and `omb cli config.json`  
- 📦 Installable with `pip` (from the repository or a single `.tar.gz` / `.whl` file)  
- 🧾 Automatic **JSON-based reproducibility** — `config.json` is in the output folder **and** in the download ZIP  
- 🔁 Multi-replica system generation  
- 📦 Self-contained outputs with input files, topologies, mdp files and run scripts  

### 🆕 New in v2.0.0

**Installation and commands**
- `pip install` support and the `omb` command (`omb gui`, `omb cli`, `omb --version`)
- Builds are written to the folder you run `omb` from, or to `$OMB_WORKDIR`

**Protein insertion page**
- **System view** of box, membrane and proteins in three projections (x–z, y–z, x–y),
  drawn from the real protein beads and rotated exactly as COBY will place them
- **Box size and shape** can be adjusted directly on this page
- Warnings when a protein sticks out of the box or when proteins overlap

**Protein names**
- Every protein gets a **required name** (e.g. `1UBQ`), used for its files —
  `1UBQ.pdb`, `1UBQ.itp` — and for its moleculetype in the topology
- Names are checked (letters, digits, `_`, `-`; no duplicates)

**Multiple copies of a protein**
- **Move together** — copies on a regular grid, sharing one position and orientation (as before)
- **Move independently** — every copy has its own position and rotation

**Cleaner topology**
- Copies of the same protein share **one** include and **one** line in `[ molecules ]`:

  ```text
  #include "../../toppar/1UBQ.itp"
  ...
  [ molecules ]
  1UBQ        3
  ```

- The protein itp uses the protein name as moleculetype (`[ moleculetype ]` → `1UBQ 1`)

**Fixes**
- `config.json` is now included in the download ZIP; `omb cli` also writes it to the output folder
- `index.sh` reads the current `config.json` format, so `run.sh` creates the index groups
  (`MEMB`, `SOLV`, `PROT`, `SOLU`) without errors
- Protein rotations in the runner follow COBY's convention (Rx·Ry·Rz)

---

## 🧭 Workflow

```text
GUI → Config → Build → config.json → ZIP → CLI reuse
```

---

## 📦 Installation

OpenMembraneBuilder requires **Python 3.10 or newer** and installs with `pip`, which gives you
the `omb` command. All Python dependencies (COBY, MDAnalysis, Streamlit, vermouth/martinize2, …)
are installed automatically.

### From the compressed file

```bash
pip install openmembranebuilder-2.0.0.tar.gz
# or the wheel
pip install openmembranebuilder-2.0.0-py3-none-any.whl
```

### From the repository

```bash
git clone -b development https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
pip install .            # or: pip install -e .   (editable, for development)
```

### With conda

```bash
conda env create -f environment.yml     # installs the dependencies and OpenMembraneBuilder
conda activate OpenMembraneBuilder
```

### Optional

- **PyMOL** — only needed for protein mutations:
  `conda install -c conda-forge pymol-open-source`
- **jq** — used by `run.sh` to create the index file:
  `conda install -c conda-forge jq` (or `brew install jq` / `apt install jq`)
- **GROMACS** — to run the simulations

### Building the compressed file yourself

```bash
pip install build
python -m build          # dist/openmembranebuilder-2.0.0.tar.gz and .whl
```

---

## 🚀 Usage

### GUI mode

```bash
omb gui                      # opens the app in the browser
omb gui --port 8502          # extra options are passed to Streamlit
```

---

### CLI mode

```bash
omb cli config.json
```

Use any `config.json` from an earlier build to rebuild the same system or, after editing it,
a variant of it.

Outputs are written to the folder you run `omb` from:

```bash
OMB_WORKDIR=/scratch/me/membranes omb cli config.json   # write them somewhere else
```

`omb --version` prints the installed version.

> The previous commands `streamlit run app.py` and `python app.py --no-gui config.json`
> still work from the repository folder.

---

## 📁 Output Structure

Each build creates a self-contained project in `output/`, and a ZIP of it in `downloads/`:

```text
output/OpenMembraneBuilder-xxxx/
├── config.json          the full build configuration
├── user_inputs/         your input files (e.g. 1UBQ.pdb, 1UBQ.itp)
├── toppar/              force field and molecule itps (e.g. 1UBQ.itp)
├── mdp/                 em, eq1, eq2, eq3, md
├── scripts/             index.sh
└── systems/
    ├── run.sh           minimisation, equilibration and production for every replica
    ├── R0001/           system.gro, system.pdb, topol.top, placed proteins
    └── R0002/ …

downloads/OpenMembraneBuilder-xxxx.zip
```

To run the simulations:

```bash
cd output/OpenMembraneBuilder-xxxx/systems
./run.sh
```

✔ Fully reproducible  
✔ Ready for simulation  
✔ ZIP archive generated automatically  

---

## 🔁 Reproducibility

Every build generates:

```text
config.json
```

It is stored in the output folder and in the ZIP, and allows:

- exact rebuilds (`omb cli config.json`)
- parameter editing  
- automated workflows  

---

## 📚 Documentation

👉 **Online docs (recommended):**  
https://lolicato.github.io/OpenMembraneBuilder  

👉 Local preview:

```bash
pip install ".[docs]"
mkdocs serve
```

---

## 🤝 Contributing

We welcome:

- Bug reports  
- Feature requests  
- Documentation improvements  
- New workflows  

---

## 📄 Licensing

    OpenMembraneBuilder is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).

    Copyright (C) 2026  Fabio Lolicato

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

---

## 🙏 Acknowledgements

OpenMembraneBuilder builds upon:

- COBY  
- MDAnalysis  
- martinize2 / vermouth  
- Streamlit and streamlit-molstar  
- NumPy / SciPy  
