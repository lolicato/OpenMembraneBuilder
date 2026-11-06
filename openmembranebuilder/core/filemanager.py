"""Copying files into a build's output folder and packing the result into a zip archive.

Every build gets the same set of supporting files next to the systems: the force field, the
simulation settings, the run scripts and copies of the files you uploaded. The `FileManager`
puts them in place and, at the end, packs the whole folder into `downloads/<folder>.zip`.
"""

from __future__ import annotations
import sys
from pathlib import Path

import shutil
import os
from openmembranebuilder.core.global_info import GlobalInfo


class FileManager():
    """Copies shipped and uploaded files into an output folder and creates the download zip."""

    def __init__(self):
        """Prepare the file manager with the standard folder locations (`GlobalInfo`)."""
        self.global_info = GlobalInfo()

    def copy_ff_folder(self, ff_name: str, destination: str):
        """Copy one force field (its folder of `.itp` files) into the output.

        An older copy in the destination is removed first, so the output always holds exactly
        the shipped version.

        Args:
            ff_name: Name of the force field, e.g. `"martini_v3"`.
            destination: The output's `toppar` folder.
        """
        src_folder = os.path.join(self.global_info.toppar_folder_path, ff_name)
        dst_folder = os.path.join(destination, ff_name)

        src_itp = os.path.join(self.global_info.toppar_folder_path, f"{ff_name}.itp")
        dst_itp = os.path.join(destination, f"{ff_name}.itp")

        if os.path.exists(dst_folder):
            shutil.rmtree(dst_folder)
        if os.path.exists(dst_itp):
            os.remove(dst_itp)

        if os.path.exists(src_folder):
            shutil.copytree(src_folder, dst_folder, dirs_exist_ok=True)
        if os.path.exists(src_itp):
            shutil.copy(src_itp, dst_itp)

    def copy_mdp_files(self, ff_name: str, destination: str, system_type: str):
        """Copy the simulation settings (`em.mdp`, `eq1.mdp`, … `md.mdp`) into `<output>/mdp`.

        Args:
            ff_name: Name of the force field, e.g. `"martini_v3"`.
            destination: The output folder.
            system_type: `"membrane"` for membrane-only systems, `"protein"` when the system
                contains proteins (these settings couple protein and membrane together).
        """
        src_mdp = os.path.join(self.global_info.mdp_folder_path, ff_name, system_type)
        dst_mdp = os.path.join(destination, "mdp")
        if os.path.exists(src_mdp):
            os.makedirs(dst_mdp, exist_ok=True)
            for file in os.listdir(src_mdp):
                shutil.copy(os.path.join(src_mdp, file), os.path.join(dst_mdp, file))

    def copy_scripts(self, destination: str):
        """Copy the run scripts into the output.

        `run.sh` goes into `<output>/systems` (it runs every replica found there) and the
        helper `index.sh`, which creates the index groups, goes into `<output>/scripts`.

        Args:
            destination: The output folder.
        """
        src_scripts = self.global_info.scripts_folder_path
        dst_scripts = os.path.join(destination, "systems")

        if os.path.exists(src_scripts):
            shutil.copytree(src_scripts, dst_scripts, dirs_exist_ok=True)

        # Move helper scripts to output/scripts
        src_index = os.path.join(dst_scripts, "index.sh")
        dst_index_dir = os.path.join(destination, "scripts")

        if os.path.exists(src_index):
            os.makedirs(dst_index_dir, exist_ok=True)
            shutil.move(src_index, os.path.join(dst_index_dir, "index.sh"))

    def create_zip_folder(self, folder_path: str) -> str:
        """Pack a whole output folder into a zip archive in the downloads folder.

        Args:
            folder_path: The output folder of the build.

        Returns:
            The full path of the zip file, e.g. `downloads/OpenMembraneBuilder-<time>.zip`.
        """
        folder_path = Path(folder_path)
        downloads_dir = Path(self.global_info.downloads_folder_path)
        downloads_dir.mkdir(parents=True, exist_ok=True)

        zip_base = downloads_dir / folder_path.name
        shutil.make_archive(str(zip_base), 'zip', str(folder_path))

        return str(zip_base.with_suffix('.zip'))

    @staticmethod
    def add_to_zip(zip_path: str, file_path: str, arcname: str | None = None) -> None:
        """Add one file to an existing zip archive, replacing an older copy of it.

        Used for `config.json`, which is written after the zip has been created.
        Nothing happens when the zip or the file does not exist.

        Args:
            zip_path: The zip archive to update.
            file_path: The file to add.
            arcname: Name of the file inside the zip; by default its own file name.
        """
        import zipfile
        if not (zip_path and os.path.exists(zip_path) and os.path.exists(file_path)):
            return
        arcname = arcname or os.path.basename(file_path)
        with zipfile.ZipFile(zip_path) as zf:
            present = arcname in zf.namelist()
        if present:  # rewrite the zip without the old copy
            tmp = zip_path + ".tmp"
            with zipfile.ZipFile(zip_path) as src, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as dst:
                for item in src.infolist():
                    if item.filename != arcname:
                        dst.writestr(item, src.read(item.filename))
            os.replace(tmp, zip_path)
        with zipfile.ZipFile(zip_path, "a", zipfile.ZIP_DEFLATED) as zf:
            zf.write(file_path, arcname)

    def copy_userinput(self, file_path, user_inputs_dir):
        """Keep a copy of a file you provided in the output's `user_inputs` folder.

        Args:
            file_path: The uploaded or prepared file.
            user_inputs_dir: The output's `user_inputs` folder (created when missing).

        Returns:
            The path of the copy.
        """
        os.makedirs(user_inputs_dir, exist_ok=True)
        new_file_path = os.path.join(user_inputs_dir, os.path.basename(file_path))
        if os.path.abspath(file_path) != os.path.abspath(new_file_path):
            shutil.copy2(file_path, new_file_path)
        return new_file_path
