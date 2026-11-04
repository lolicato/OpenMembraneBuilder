import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import shutil
import os
from src.core.global_info import GlobalInfo


class FileManager():
    def __init__(self):
        self.global_info = GlobalInfo()

    def copy_ff_folder(self, ff_name: str, destination: str):
        '''copies all the forcefield files into the output folder'''
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
        '''copies mdp files into the output folder'''
        src_mdp = os.path.join(self.global_info.mdp_folder_path, ff_name, system_type)
        dst_mdp = os.path.join(destination, "mdp")
        if os.path.exists(src_mdp):
            os.makedirs(dst_mdp, exist_ok=True)
            for file in os.listdir(src_mdp):
                shutil.copy(os.path.join(src_mdp, file), os.path.join(dst_mdp, file))

    def copy_scripts(self, destination: str):
        """Copy scripts into the output folder."""
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
        '''creates a downloadable zip folder from the output (absolute path, independent of cwd)'''
        folder_path = Path(folder_path)
        downloads_dir = Path(self.global_info.downloads_folder_path)
        downloads_dir.mkdir(parents=True, exist_ok=True)

        zip_base = downloads_dir / folder_path.name
        shutil.make_archive(str(zip_base), 'zip', str(folder_path))

        return str(zip_base.with_suffix('.zip'))

    def copy_userinput(self, file_path, user_inputs_dir):
        '''copies a given user input to its destination and returns the new path'''
        os.makedirs(user_inputs_dir, exist_ok=True)
        new_file_path = os.path.join(user_inputs_dir, os.path.basename(file_path))
        if os.path.abspath(file_path) != os.path.abspath(new_file_path):
            shutil.copy2(file_path, new_file_path)
        return new_file_path
