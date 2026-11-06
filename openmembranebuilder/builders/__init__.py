"""The two builders the wizard chains together.

- `martinize` – prepares the proteins: takes coarse-grained files you already have, or turns
  an atomistic structure (uploaded or downloaded from the RCSB PDB) into a Martini model with
  martinize2.
- `membrane` – builds the systems with COBY: membrane composition, box, water and ions,
  protein placement, replicas, topology and the files to run GROMACS.

Each builder has the same four parts: `config` (its settings), `gui` (the widgets you see),
`page` (what happens when you click) and `runner` (the actual work, also used by `omb cli`).
"""

