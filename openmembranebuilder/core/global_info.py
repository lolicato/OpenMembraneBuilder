"""Where OpenMembraneBuilder finds its own files and where it writes yours.

Two kinds of locations are used:

- **Shipped data** – force fields (`toppar`), simulation settings (`mdp`), run scripts,
  membrane templates and example files. They live inside the installed package
  (`openmembranebuilder/resources/`), so they are found wherever the package is installed.
- **Your results** – the `output/` folder with the built systems, the `downloads/` folder with
  the zip archives, and a temporary folder for uploaded files. They are created in the working
  folder: the folder where `omb gui` / `omb cli` was started, or the folder named by the
  environment variable `OMB_WORKDIR`.
"""

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
    """The fixed folder and file locations used throughout the program.

    Every part of OpenMembraneBuilder creates its own `GlobalInfo()` and reads the paths it
    needs from it, so changing a location here changes it everywhere.

    Attributes:
        version: Version label of the program, used in titles.
        resources_path: Folder with all shipped data (`openmembranebuilder/resources`).
        toppar_folder: Name of the topology folder inside each output (`toppar`).
        base_folder: Folder in which every build creates its own sub-folder (`./output`).
        mdp_folder: Name of the simulation-settings folder inside each output (`mdp`).
        scripts_folder: Name of the scripts folder inside each output (`scripts`).
        membrane_template_example_folder: Name of the folder with the example template.
        membrane_templates_path: Folder with the predefined membrane templates.
        toppar_folder_path: Folder with the shipped force fields (`.top` files and their
            `.itp` folders).
        mdp_folder_path: Folder with the shipped simulation settings, per force field and
            system type (`membrane` / `protein`).
        scripts_folder_path: Folder with the shipped `run.sh` and `index.sh`.
        membrane_template_example_path: Folder with `example_template.ob`, offered for
            download on the membrane page.
        downloads_folder_path: Folder where the zip archive of every build is written.
        temp_folder: Folder for uploads during a session (removed after the build).
        default_template_name: Name used for the composition when no template with named
            blocks is used (`single_setup`).
        import_folder: Name of the folder with extra structures for Martini 2.
        chol_file: Cholesterol structure used with Martini 2.2 (`CHOL`).
        chol2_file: Improved cholesterol structure used with Martini 2.2 (`CHOL2`).
        wf_file: Antifreeze water structure used with Martini 2.2 (`WF`).
    """
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
