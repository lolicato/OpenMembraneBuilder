"""Rebuilding systems from a saved `config.json`, without the graphical interface.

This is what `omb cli config.json` runs. The configuration lists the builders in the order the
wizard used them (`martinize`, then `membrane`); each one is created from its own section of
the file and run.

Proteins and custom lipids were already prepared when the configuration was first written:
their finished files are recorded in the `membrane` section. The builders that produced them
are therefore skipped, and the membrane builder uses the recorded files directly.
"""

import importlib
import inspect
import json
import os
from datetime import datetime


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


def run_from_config(config_path: str):
    """Build the systems described in a `config.json`.

    A new output folder `output/OpenMembraneBuilder-<date>_<time>-headless` is created
    (unless the file names one), the systems are built in it, the configuration is copied
    next to them and added to the zip archive in `downloads/`.

    Args:
        config_path: Path of the `config.json`, usually one written by an earlier build.
    """
    with open(config_path) as f:
        data = json.load(f)

    # The saved membrane config already contains the RESULTS of the earlier builders
    # (protein_files = martinize outputs, molecule_builder_index = cglipid outputs),
    # and the temp folders those builders used no longer exist. Re-running them would
    # fail (or just duplicate work), so they are skipped.
    membrane = data.get("membrane", {})
    skip = set()
    if membrane.get("protein_files"):
        skip.add("martinize")
    if membrane.get("molecule_builder_index"):
        skip.add("cglipid")

    for builder_name, config_dict in data.items():
        if builder_name in skip:
            print(f"[headless] skipping '{builder_name}': its results are stored in the membrane config")
            continue

        config_module = importlib.import_module(f"openmembranebuilder.builders.{builder_name}.config")
        config = config_module.Config(**config_dict)

        if builder_name == "membrane" and not config.base_folder:
            from openmembranebuilder.core.global_info import GlobalInfo
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            config.base_folder = os.path.join(GlobalInfo().base_folder, f"OpenMembraneBuilder-{ts}-headless")

        runner_class = _load_runner_class(builder_name)
        result = runner_class().run(config)

        # keep the configuration with the output and in the download zip, as the GUI does
        if builder_name == "membrane" and getattr(config, "base_folder", None):
            from openmembranebuilder.core.filemanager import FileManager
            out_cfg = os.path.join(config.base_folder, "config.json")
            if os.path.abspath(out_cfg) != os.path.abspath(config_path):
                with open(out_cfg, "w") as f:
                    json.dump(data, f, indent=4)
            if isinstance(result, tuple) and len(result) == 2:
                FileManager.add_to_zip(result[1], out_cfg, "config.json")
