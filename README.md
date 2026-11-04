# 🚀 OpenMembraneBuilder-v2.0.0

![Version](https://img.shields.io/badge/version-v2.0.0-green)
![Python](https://img.shields.io/badge/python-3.10-blue)
![Status](https://img.shields.io/badge/status-active-success)

👉 Working Progress!


🌐 **📖 Documentation:** https://lolicato.github.io/OpenMembraneBuilder  
👉 *Full user guide, and API reference available online*

---

**OpenMembraneBuilder** is a modern **GUI + CLI toolkit for building MARTINI coarse-grained membrane and membrane–protein systems**, powered by COBY.

It combines:
- 🧬 **MDAnalysis** for structure handling  
- ⚙️ **COBY** for system generation  
- 🌐 **Streamlit** for an interactive GUI  

with a strong focus on **reproducibility**, **usability**, and **local execution**.

---

## 🖼️ Interface

![Application GUI](./pictures/GUI.png)

---

## ✨ Features

- 🧪 Build membrane and membrane–protein systems  
- 🎛️ Fully interactive Streamlit GUI  
- 🧾 Automatic **JSON-based reproducibility**  
- 💻 Command-line interface (CLI) for batch workflows  
- 📦 Self-contained outputs with input files included  
- 🔁 Multi-replica system generation  

---

## 🧭 Workflow

```text
GUI → Config → Build → config.json → ZIP → CLI reuse
```

---

## 📦 Installation

### Clone the Development version


```bash
git clone https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
```

### Create environment

```bash
conda env create -f environment.yml
conda activate OpenMembraneBuilder
```

---

## 🚀 Usage

### GUI mode

```bash
streamlit run app.py
```

---

### CLI mode

```bash
python app.py --no-gui config.json
```

---

## 📁 Output Structure

Each build creates a self-contained project:

```text
outputs/OpenMembraneBuilder-xxxx/
├── config.json
├── user_inputs/
├── simulations/
├── toppar/
└── mdp/
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

This allows:

- exact rebuilds  
- parameter editing  
- automated workflows  

---

## 📚 Documentation

👉 **Online docs (recommended):**  
https://lolicato.github.io/OpenMembraneBuilder  

👉 Local preview:

```bash
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
- Streamlit  
- NumPy / SciPy  

