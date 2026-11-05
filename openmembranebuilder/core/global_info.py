from dataclasses import dataclass
import os
from pathlib import Path

# Shipped data (force fields, mdp, scripts, templates) lives inside the package, so it is found
# both in a source checkout and after `pip install`.
_PKG = Path(__file__).resolve().parents[1]          # .../openmembranebuilder
_RESOURCES = str(_PKG / "resources")
# Everything the program writes (outputs, zip downloads, uploads) goes to the working folder:
# the folder where `omb gui` / `omb cli` is started, or $OMB_WORKDIR.
_WORK = Path(os.environ.get("OMB_WORKDIR") or Path.cwd()).resolve()


@dataclass
class GlobalInfo:
    """Absolute paths: shipped data from the package, outputs in the working folder."""
    version: str = "v1.2.1"
    resources_path: str = _RESOURCES
    toppar_folder: str = "toppar"
    base_folder: str = str(_WORK / "output")
    mdp_folder: str = "mdp"
    scripts_folder: str = "scripts"
    membrane_template_example_folder: str = "example_files"
    membrane_templates_path: str = os.path.join(_RESOURCES, "membrane_templates")
    toppar_folder_path: str = os.path.join(_RESOURCES, "toppar")
    mdp_folder_path: str = os.path.join(_RESOURCES, "mdp")
    scripts_folder_path: str = os.path.join(_RESOURCES, "scripts")
    membrane_template_example_path: str = os.path.join(_RESOURCES, "example_files")
    downloads_folder_path: str = str(_WORK / "downloads")
    temp_folder: str = str(_WORK / "temp_uploads")
    default_template_name: str = "single_setup"
    import_folder: str = "file_imports"
    chol_file: str = os.path.join(_RESOURCES, "file_imports", "martini_v2.2", "CHOL.pdb")
    chol2_file: str = os.path.join(_RESOURCES, "file_imports", "martini_v2.2", "CHOL2.pdb")
    wf_file: str = os.path.join(_RESOURCES, "file_imports", "martini_v2.2", "WF.gro")
