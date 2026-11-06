"""Building systems from a `config.json`, without the graphical interface (`omb cli`).

The configuration has one section per builder, in the order the wizard used them:

- `martinize` – how the proteins are prepared (download or file, clean-up, chains, mutation,
  martinize2 options);
- `membrane` – the system: box, water, lipids, proteins and their placement.

**Proteins.** A configuration written by an earlier build already lists the finished protein
files (`membrane.protein_files`); they are used as they are, and martinize2 is not run again.
The proteins are prepared (again) when

- the configuration has a `martinize` section but no finished protein files,
- `--martinize` is given, or
- `--pdb` is given: it replaces the protein(s) of the configuration by other structures – a
  PDB code to download from RCSB, or an atomistic `.pdb` file – which are then prepared with
  the martinize2 settings of the configuration. This turns any configuration into a template
  for many proteins:

  ```bash
  for PDB in 1BFF 1UBQ 2RH1; do omb cli template.json --pdb $PDB; done
  ```

The work of the protein preparation is kept in `<output>/martinize/`, and the configuration
written to the output records the prepared files, so the output rebuilds without running
martinize2 again.
"""

import copy
import importlib
import inspect
import json
import os
import re
import shutil
from datetime import datetime

_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def _load_runner_class(builder_name: str):
    """Find the builder class (`...Runner`) of one builder.

    Args:
        builder_name: Name of the builder, e.g. `"membrane"` or `"martinize"`.

    Returns:
        The runner class defined in `openmembranebuilder.builders.<builder_name>.runner`.

    Raises:
        ImportError: The builder has no runner class.
    """
    module = importlib.import_module(f"openmembranebuilder.builders.{builder_name}.runner")
    for _, obj in vars(module).items():
        # only classes DEFINED in this module (not imported helpers)
        if (
            inspect.isclass(obj)
            and obj.__module__ == module.__name__
            and obj.__name__.endswith("Runner")
        ):
            return obj
    raise ImportError(f"No *Runner class found in openmembranebuilder.builders.{builder_name}.runner")


def _apply_pdb_overrides(martinize: dict, pdbs: list) -> list:
    """Put other structures into the protein settings of a configuration.

    The n-th `--pdb` replaces the source of the n-th protein; its other settings
    (clean-up, martinize2 force field, elastic network, …) stay as they are. Extra `--pdb`
    values get a copy of the last protein's settings. Chain and residue selections and
    mutations belong to one particular structure, so they are reset.

    Args:
        martinize: The `martinize` section of the configuration, changed in place.
        pdbs: PDB codes (downloaded from RCSB) or paths of atomistic `.pdb` files.

    Returns:
        The protein names, in order (the PDB code, or the file name without `.pdb`).

    Raises:
        ValueError: A name that cannot be used as a protein name.
    """
    proteins = martinize.setdefault("proteins", [])
    if not proteins:
        proteins.append({"forcefield": "martini3001", "clean_structure": True})
    names = []
    for i, src in enumerate(pdbs):
        if i >= len(proteins):
            proteins.append(copy.deepcopy(proteins[-1]))
        p = proteins[i]
        if os.path.isfile(src):                                   # an atomistic file
            name = os.path.splitext(os.path.basename(src))[0]
            p.update(atomistic_source="Upload PDB", atomistic_pdb_path=os.path.abspath(src),
                     fetch_pdb_id="")
        else:                                                     # a PDB code
            name = src.strip().upper()
            p.update(atomistic_source="Fetch from PDB database", fetch_pdb_id=name,
                     atomistic_pdb_path="")
        if not _NAME_RE.match(name):
            raise ValueError(f"'{name}' cannot be used as a protein name "
                             "(letters, digits, _ and - only)")
        p.update(name=name, protein_input_mode="martinize",
                 selected_chains=[], chain_ranges={}, do_mutation=False,
                 output_cg_pdb_path="", output_itp_path="", base_folder="")
        names.append(name)
    martinize["build_mode"] = "protein_membrane"
    martinize["n_proteins"] = len(proteins)
    return names


def _prepare_proteins(martinize: dict, work_dir: str):
    """Prepare every protein of a `martinize` section, as the protein-preparation pages do.

    Args:
        martinize: The `martinize` section of the configuration.
        work_dir: Folder for the work of each protein (`<work_dir>/p0`, `p1`, …).

    Returns:
        `(config, protein_files)`: the protein settings with their finished files, and the
        proteins as the membrane builder expects them (`{"1UBQ.pdb": {...}}`).

    Raises:
        RuntimeError: A protein could not be prepared (the message says why).
    """
    from openmembranebuilder.builders.martinize.config import Config as MartinizeConfig
    from openmembranebuilder.builders.martinize.runner import MartinizeRunner

    mconfig = MartinizeConfig(**martinize)
    mconfig.n_proteins = len(mconfig.proteins)
    runner = MartinizeRunner()
    for i, pc in enumerate(mconfig.proteins):
        label = pc.name or f"protein {i + 1}"
        if not pc.name or not _NAME_RE.match(pc.name):
            raise RuntimeError(f"{label}: every protein needs a valid 'name' "
                               "(letters, digits, _ and - only)")
        pc.base_folder = os.path.join(work_dir, f"p{i}")
        print(f"[headless] preparing {label} ...", flush=True)
        result = runner.run(mconfig, protein_index=i) or {}
        if not result.get("success"):
            raise RuntimeError(f"{label}: {result.get('message') or 'preparation failed'}")
        pc.output_cg_pdb_path = result["outputs"].get("cg_pdb_path", "")
        pc.output_itp_path = result["outputs"].get("itp_path", "")
        print(f"[headless] {label}: {pc.output_cg_pdb_path}", flush=True)
    return mconfig, mconfig.protein_files()


def _rekey_poses(membrane: dict, old_names: list, new_names: list) -> None:
    """Keep the protein placement of the configuration when the proteins get new names.

    Positions and rotations are stored per protein file name (`"1UBQ.pdb"`). When the
    proteins are replaced (`--pdb`), the n-th old protein's placement is handed to the n-th
    new protein; placements of proteins that no longer exist are dropped.

    Args:
        membrane: The `membrane` section, changed in place.
        old_names: Protein file names before, in order.
        new_names: Protein file names after, in order.
    """
    mapping = dict(zip(old_names, new_names))
    for key in ("protein_params",):
        per_system = membrane.get(key) or {}
        for system, poses in per_system.items():
            per_system[system] = {mapping[k]: v for k, v in poses.items() if k in mapping}
    ins = membrane.get("insertion_params") or {}
    membrane["insertion_params"] = {mapping[k]: v for k, v in ins.items() if k in mapping}


def run_from_config(config_path: str, pdbs: list = None, martinize: bool = False):
    """Build the systems described in a `config.json`.

    A new output folder `output/OpenMembraneBuilder-<date>_<time>-headless` is created (with
    `--pdb`, the protein names are added to it), the proteins are prepared when needed (see
    the module description), the systems are built, and the configuration – with the
    prepared protein files and the placements actually used – is written next to them and
    added to the zip archive in `downloads/`.

    Args:
        config_path: Path of the `config.json` (written by an earlier build, or a template).
        pdbs: Other structures for the proteins: PDB codes or atomistic `.pdb` files.
        martinize: Prepare the proteins again even when finished files are recorded.

    Raises:
        RuntimeError: A protein could not be prepared.
    """
    with open(config_path) as f:
        data = json.load(f)

    from openmembranebuilder.core.global_info import GlobalInfo
    from openmembranebuilder.core.filemanager import FileManager

    membrane = data.setdefault("membrane", {})
    old_names = list((membrane.get("protein_files") or {}).keys())

    names = []
    if pdbs:
        names = _apply_pdb_overrides(data.setdefault("martinize", {}), pdbs)
        martinize = True

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = "headless" + ("-" + "-".join(names) if names else "")
    base_folder = os.path.join(GlobalInfo().base_folder, f"OpenMembraneBuilder-{ts}-{tag}")

    # --- proteins: use the recorded files, or prepare them -----------------------------------
    mdata = data.get("martinize") or {}
    if mdata.get("proteins") and (martinize or not membrane.get("protein_files")):
        try:
            mconfig, protein_files = _prepare_proteins(mdata, os.path.join(base_folder, "martinize"))
        except RuntimeError as e:
            # keep the work folder only when it holds something to look at (e.g. martinize.log)
            has_files = any(files for _, _, files in os.walk(base_folder))
            if has_files:
                raise RuntimeError(f"{e}\n  work files: {os.path.join(base_folder, 'martinize')}") from None
            shutil.rmtree(base_folder, ignore_errors=True)
            raise
        data["martinize"] = mconfig.to_dict()
        _rekey_poses(membrane, old_names, list(protein_files.keys()))
        membrane["protein_files"] = protein_files
    elif membrane.get("protein_files") and "martinize" in data:
        print("[headless] using the protein files recorded in the configuration "
              "(--martinize prepares them again)")

    # --- other builders (e.g. custom lipids): their results are recorded in 'membrane' ------
    for builder_name, config_dict in data.items():
        if builder_name in ("martinize", "membrane"):
            continue
        if builder_name == "cglipid" and membrane.get("molecule_builder_index"):
            print(f"[headless] skipping '{builder_name}': its results are stored in the membrane config")
            continue
        config_module = importlib.import_module(f"openmembranebuilder.builders.{builder_name}.config")
        _load_runner_class(builder_name)().run(config_module.Config(**config_dict))

    # --- the systems --------------------------------------------------------------------------
    config_module = importlib.import_module("openmembranebuilder.builders.membrane.config")
    config = config_module.Config(**membrane)
    if not config.base_folder:
        config.base_folder = base_folder
    result = _load_runner_class("membrane")().run(config)

    # keep the configuration with the output and in the download zip, as the GUI does;
    # it records the prepared proteins and the placements that were used
    membrane["protein_files"] = config.protein_files
    membrane["protein_params"] = config.protein_params
    out_cfg = os.path.join(config.base_folder, "config.json")
    if os.path.abspath(out_cfg) != os.path.abspath(config_path):
        with open(out_cfg, "w") as f:
            json.dump(data, f, indent=4, default=str)
    if isinstance(result, tuple) and len(result) == 2:
        FileManager.add_to_zip(result[1], out_cfg, "config.json")
    print(f"[headless] done: {config.base_folder}")
