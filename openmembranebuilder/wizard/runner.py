import importlib
import inspect
import json
import os
from datetime import datetime


def _load_runner_class(builder_name: str):
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
