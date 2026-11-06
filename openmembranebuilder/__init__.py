"""OpenMembraneBuilder: build Martini coarse-grained membranes and membrane–protein systems.

The package is organised in four parts:

- `openmembranebuilder.cli` – the `omb` command (`omb gui`, `omb cli config.json`).
- `openmembranebuilder.wizard` – the step-by-step graphical workflow (welcome page,
  protein preparation, membrane, results) and the headless rebuild from a `config.json`.
- `openmembranebuilder.builders` – the two builders the wizard chains together:
  `martinize` prepares coarse-grained proteins, `membrane` builds the systems with COBY.
- `openmembranebuilder.core` – shared helpers: folder locations, file copying, page layout.

Force fields, simulation settings (mdp files), run scripts and membrane templates are shipped
inside the package, in `openmembranebuilder/resources/`.
"""

