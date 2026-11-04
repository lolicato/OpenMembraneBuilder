from dataclasses import dataclass
import os
from pathlib import Path

# src/core/global_info.py -> project root
_ROOT = Path(__file__).resolve().parents[2]
_RESOURCES = str(_ROOT / "resources")


@dataclass
class GlobalInfo:
    """All paths are absolute now, so the app no longer depends on the working directory."""
    version: str = "v1.2.1"
    resources_path: str = _RESOURCES
    toppar_folder: str = "toppar"
    base_folder: str = str(_ROOT / "output")
    mdp_folder: str = "mdp"
    scripts_folder: str = "scripts"
    membrane_template_example_folder: str = "example_files"
    membrane_templates_path: str = os.path.join(_RESOURCES, "membrane_templates")
    toppar_folder_path: str = os.path.join(_RESOURCES, "toppar")
    mdp_folder_path: str = os.path.join(_RESOURCES, "mdp")
    scripts_folder_path: str = os.path.join(_RESOURCES, "scripts")
    membrane_template_example_path: str = os.path.join(_RESOURCES, "example_files")
    downloads_folder_path: str = str(_ROOT / "downloads")
    temp_folder: str = str(_ROOT / "temp_uploads")
    default_template_name: str = "single_setup"
    import_folder: str = "file_imports"
    chol_file: str = os.path.join(_RESOURCES, "file_imports", "martini_v2.2", "CHOL.pdb")
    chol2_file: str = os.path.join(_RESOURCES, "file_imports", "martini_v2.2", "CHOL2.pdb")
    wf_file: str = os.path.join(_RESOURCES, "file_imports", "martini_v2.2", "WF.gro")
