import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[3]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import os
import re
import math
import shutil
import random
import numpy as np
import MDAnalysis as mda
from datetime import datetime
from scipy.spatial.transform import Rotation as R
from typing import Optional
import COBY


from src.core.global_info import GlobalInfo
from src.core.filemanager import FileManager
from src.builders.membrane.parser import MembraneParser



class MembraneRunner:

    def __init__(self):
        self.global_info      = GlobalInfo()
        self.file_manager     = FileManager()
        self.topology_editor  = TopologyEditor()
        self.membrane_parser  = MembraneParser()

    # ─── Public entry point ───────────────────────────────────────────────────

    def run(self, config) -> list[str]:
        """
        Build all membrane systems defined in config.
        Handles both membrane-only and protein insertion.
        Returns (output_pdbs, zip_path).
        """
        with_protein = bool(config.protein_files)

        base_folder    = config.base_folder or config.output_name
        systems_folder = os.path.join(base_folder, "systems")
        toppar_folder  = os.path.join(base_folder, self.global_info.toppar_folder)
        user_inputs_dir  = os.path.join(base_folder, "user_inputs")
        os.makedirs(base_folder,    exist_ok=True)
        os.makedirs(systems_folder, exist_ok=True)
        os.makedirs(toppar_folder,  exist_ok=True)
        os.makedirs(user_inputs_dir,   exist_ok=True)

        # Step 1 — copy shared files
        mdp_target = "protein" if with_protein else "membrane"
        self.file_manager.copy_ff_folder(config.selected_ff, toppar_folder)
        self.file_manager.copy_mdp_files(config.selected_ff, base_folder, mdp_target)
        self.file_manager.copy_scripts(base_folder)

        molecule_imports = self._build_molecule_imports(config)
        itp_input_ff     = f"include:{self.global_info.toppar_folder_path}/{config.selected_ff}.top"

        if with_protein:
            for pdb_name, paths in config.protein_files.items():
                original_itp = paths["itp_path"]
                stem = os.path.splitext(pdb_name)[0]
                itp_dest = os.path.join(toppar_folder, f"{stem}.itp")
                if os.path.abspath(original_itp) != os.path.abspath(itp_dest):
                    shutil.copy2(original_itp, itp_dest)
                moleculetype_name = stem if len(config.protein_files) > 1 else "Protein"
                self.topology_editor.overwrite_moleculetype_line(itp_dest, moleculetype_name=moleculetype_name)
                paths["itp_path"] = itp_dest
                self.file_manager.copy_userinput(original_itp, user_inputs_dir)
        if config.template_active and config.template_path and os.path.exists(config.template_path):              # ← new block
                config.template_path = self.file_manager.copy_userinput(config.template_path, user_inputs_dir)
        # Templates — entries is either a dict {name: entries} or a flat list
        templates = config.entries if isinstance(config.entries, dict) else {
            self.global_info.default_template_name: config.entries
        }
        if with_protein:
            for pdb_name, paths in config.protein_files.items():
                src = paths["pdb_path"]
                if os.path.exists(src):
                    new_path = os.path.join(user_inputs_dir, pdb_name)
                    if os.path.abspath(src) != os.path.abspath(new_path):
                        paths["pdb_path"] = self.file_manager.copy_userinput(src, user_inputs_dir)

        all_systems = []

        for template_name, entries in templates.items():
            
            config.entries_current = entries                                          
            config.membrane_string = self.membrane_parser.create_membrane_str(config) if entries else ""
            selected_molecule_builders = self._collect_selected_molecule_builders(config)
            suffix = f"{template_name}_" if template_name != self.global_info.default_template_name else ""

            params = {
                "boxx":               float(config.box_x),
                "boxy":               float(config.box_y),
                "boxz":               float(config.box_z),
                "box_type":           config.box_type,
                "moleculeimport":     molecule_imports,
                "molecule_builder":   selected_molecule_builders,
                "membrane":           config.membrane_string,
                "solvation":          config.solvation,
                "selectedforcefield": config.selected_ff,
            }

            # Step 2 — per-system pipeline
            for i in range(config.n_systems):
                folder_name = f"{suffix}R{i + 1:04d}"
                folder      = os.path.join(systems_folder, folder_name)
                os.makedirs(folder, exist_ok=True)

                if with_protein:
                    itp_includes = [f"include:{paths['itp_path']}" for paths in config.protein_files.values()]

                    if config.more_copies and config.copy_number > 1:
                        params["itp_input"] = [itp_input_ff] + itp_includes * config.copy_number
                        # _insert_multiple_proteins still works on a single protein — see note below
                        first_pdb_name = next(iter(config.protein_files))
                        first_pdb_dest = os.path.join(folder, first_pdb_name)
                        shutil.copy2(config.protein_files[first_pdb_name]["pdb_path"], first_pdb_dest)
                        protein_line = self._insert_multiple_proteins(
                            first_pdb_dest, folder, config, molecule_imports, system_index=i
                        )
                    else:
                        params["itp_input"] = [itp_input_ff] + itp_includes
                        protein_line = self._insert_protein(
                            folder, config, molecule_imports, system_index=i
                        )
                else:
                    protein_line        = []
                    params["itp_input"] = [itp_input_ff]

                coby_output = run_coby_simulation(params, protein_line, folder)
                self.topology_editor.edit_topology(config.selected_ff, folder)
                self.topology_editor.replace_placeholder_title(os.path.join(folder, "topol.top"), config)

                if with_protein:
                    itp_names = [
                        f"{os.path.splitext(pdb_name)[0]}.itp"
                        for pdb_name in config.protein_files
                    ]
                    self.topology_editor.fix_protein_includes_only(
                        os.path.join(coby_output, "topol.top"), itp_names
                    )
                    if config.more_copies and config.copy_number > 1:
                        self.topology_editor.fix_multiple_protein_topology(
                            os.path.join(coby_output, "topol.top")
                        )

                all_systems.append(folder)

        # Step 3 — collect PDBs for visualization
        output_pdbs = self._collect_output_pdbs(all_systems)

        zip_path = self.file_manager.create_zip_folder(base_folder)

        return output_pdbs, zip_path


    def _collect_selected_molecule_builders(self, config) -> list[str]:
        index = getattr(config, "molecule_builder_index", {}) or {}
        entries = getattr(config, "entries_current", []) or []

        selected_lipids = set()
        for entry in entries:
            if entry and len(entry) > 0:
                selected_lipids.add(entry[0])

        return [index[name] for name in sorted(selected_lipids) if name in index]
    # ─── Protein insertion ────────────────────────────────────────────────────


    def _rotate_protein(self, config, system_index: int) -> None:
        if system_index != 0 and not config.randomize_rot_every:
            return

        for pdb_name in config.protein_files:
            source = config.protein_params.get("R0001", {}).get(pdb_name, {})

            if config.randomize_rot and not config.randomize_rot_every:
                rx, ry, rz = [random.uniform(-180.0, 180.0) for _ in range(3)]
            elif not config.randomize_rot:
                rx = source.get("rx", 0.0)
                ry = source.get("ry", 0.0)
                rz = source.get("rz", 0.0)

            for i in range(config.n_systems):
                key = f"R{i + 1:04d}"
                if config.randomize_rot and config.randomize_rot_every:
                    rx, ry, rz = [random.uniform(-180.0, 180.0) for _ in range(3)]
                config.protein_params[key][pdb_name]["rx"] = round(rx, 4)
                config.protein_params[key][pdb_name]["ry"] = round(ry, 4)
                config.protein_params[key][pdb_name]["rz"] = round(rz, 4)

    def _position_protein(self, config, system_index: int) -> None:
        margin = 2.5
        if system_index != 0 and not config.randomize_pos_every:
            return

        for pdb_name in config.protein_files:
            source = config.protein_params.get("R0001", {}).get(pdb_name, {})

            if config.randomize_pos and not config.randomize_pos_every:
                cx = random.uniform(-(config.box_x / 2 - margin), config.box_x / 2 - margin)
                cy = random.uniform(-(config.box_y / 2 - margin), config.box_y / 2 - margin)
            elif not config.randomize_pos:
                cx = source.get("cx", 0.0)
                cy = source.get("cy", 0.0)

            for i in range(config.n_systems):
                key = f"R{i + 1:04d}"
                if config.randomize_pos and config.randomize_pos_every:
                    cx = random.uniform(-(config.box_x / 2 - margin), config.box_x / 2 - margin)
                    cy = random.uniform(-(config.box_y / 2 - margin), config.box_y / 2 - margin)
                config.protein_params[key][pdb_name]["cx"] = round(cx, 4)
                config.protein_params[key][pdb_name]["cy"] = round(cy, 4)

    def _resolve_cz(self, pdb_path: str, system_path: str, config,
            molecule_imports: list, rx: float, ry: float, rz: float,
            pdb_name: str = None) -> float:
        if config.z_method != "Height above Membrane":
            filename = pdb_name or os.path.basename(pdb_path)
            return config.protein_params.get("R0001", {}).get(filename, {}).get("cz", 0.0)
            
        upper_z_mem = self._create_pure_membrane_system(system_path, config, molecule_imports)
        cz = upper_z_mem / 10.0 + config.distance_to_mem - config.box_z / 2

        u   = mda.Universe(pdb_path)
        com = u.atoms.center_of_geometry()
        rot = R.from_euler("xyz", [rx, ry, rz], degrees=True)
        rotated = rot.apply(u.atoms.positions - com) + com
        cz += (rotated[:, 2].mean() - rotated[:, 2].min()) / 10.0
        return round(cz, 4)
    
    def _build_protein_line(self, pdb_path, cx, cy, cz, rx, ry, rz,
                        moleculetype: str = "Protein") -> str:
        return (
            f"file:{pdb_path} cx:{cx:.4f} cy:{cy:.4f} cz:{cz:.4f}"
            f" rx:{rx:.4f} ry:{ry:.4f} rz:{rz:.4f}"
            f" moleculetype:{moleculetype}"
        )
    
    def _insert_multiple_proteins(self, pdb_path: str, system_path: str, config,
                               molecule_imports: list, system_index: int) -> list[str]:
        self._ensure_protein_params(config)
        self._rotate_protein(config, system_index)

        pdb_name  = next(iter(config.protein_files))   # single protein for multi-copy
        key       = f"R{system_index + 1:04d}"
        system_rx = config.protein_params[key][pdb_name]["rx"]
        system_ry = config.protein_params[key][pdb_name]["ry"]
        system_rz = config.protein_params[key][pdb_name]["rz"]

        grid   = self._make_grid(config.box_x, config.box_y, config.copy_number)
        random.shuffle(grid)
        points = grid[:config.copy_number]

        cz = self._resolve_cz(pdb_path, system_path, config, molecule_imports,
                            system_rx, system_ry, system_rz)

        lines = []
        for cx, cy in points:
            if config.randomize_rot_all_copies:
                rx = round(random.uniform(-180.0, 180.0), 4)
                ry = round(random.uniform(-180.0, 180.0), 4)
                rz = round(random.uniform(-180.0, 180.0), 4)
            else:
                rx, ry, rz = system_rx, system_ry, system_rz

            lines.append(self._build_protein_line(pdb_path, cx, cy, cz, rx, ry, rz))

        return lines  # list of strings, one per copy

    def _insert_protein(self, system_path: str, config,
                    molecule_imports: list, system_index: int = 0) -> list[str]:
        self._ensure_protein_params(config)
        self._rotate_protein(config, system_index)

        # Distance mode overrides cx/cy
        if getattr(config, "use_protein_distance", False) and len(config.protein_files) == 2:
            self._position_two_proteins_by_distance(config, system_index)
        else:
            self._position_protein(config, system_index)

        key   = f"R{system_index + 1:04d}"
        lines = []

        for pdb_name, paths in config.protein_files.items():
            pdb_path = os.path.join(system_path, pdb_name)
            if not os.path.exists(pdb_path):
                shutil.copy2(paths["pdb_path"], pdb_path)

            params = config.protein_params[key][pdb_name]
            cz = self._resolve_cz(
                pdb_path, system_path, config, molecule_imports,
                params["rx"], params["ry"], params["rz"],
                pdb_name=pdb_name
            )
            params["cz"] = cz
            moleculetype = os.path.splitext(pdb_name)[0] if len(config.protein_files) > 1 else "Protein"
            lines.append(self._build_protein_line(
                pdb_path, params["cx"], params["cy"], cz,
                params["rx"], params["ry"], params["rz"],
                moleculetype=moleculetype,
            ))

        return lines
    
    def _get_protein_extent_on_x(self, pdb_path: str, rx: float, ry: float, rz: float) -> tuple[float, float]:
        """
        Returns (half_extent, center_offset) in nm.
        Rotates the protein, then returns:
        - half the x-span (max_x - min_x) / 2
        - the x offset of the COM after rotation (to handle asymmetric proteins)
        """
        u   = mda.Universe(pdb_path)
        pos = u.atoms.positions  
        com = pos.mean(axis=0)
        rot = R.from_euler("xyz", [rx, ry, rz], degrees=True)
        rotated = rot.apply(pos - com)  # centered at origin

        x_min = rotated[:, 0].min() / 10.0  
        x_max = rotated[:, 0].max() / 10.0
        half_extent = (x_max - x_min) / 2.0
        return half_extent  


    def _position_two_proteins_by_distance(self, config, system_index: int) -> None:
        """
        Positions exactly two proteins along the X axis with `config.protein_distance`
        gap between their surfaces. Both get cy=0.0.

        Protein A placed at  +x, Protein B at -x.
        x = (half_extent_A + distance/2)  and  -(half_extent_B + distance/2)
        """
        pdb_names = list(config.protein_files.keys())
        assert len(pdb_names) == 2, "Distance mode requires exactly 2 proteins"

        name_a, name_b = pdb_names

        for i in range(config.n_systems):
            key = f"R{i + 1:04d}"

            params_a = config.protein_params[key][name_a]
            params_b = config.protein_params[key][name_b]

            ext_a = self._get_protein_extent_on_x(
                config.protein_files[name_a]["pdb_path"],
                params_a["rx"], params_a["ry"], params_a["rz"],
            )
            ext_b = self._get_protein_extent_on_x(
                config.protein_files[name_b]["pdb_path"],
                params_b["rx"], params_b["ry"], params_b["rz"],
            )

            half_gap = config.protein_distance / 2.0

            x_a = +(ext_a + half_gap)
            x_b = -(ext_b + half_gap)

            config.protein_params[key][name_a]["cx"] = round(x_a, 4)
            config.protein_params[key][name_a]["cy"] = 0.0
            config.protein_params[key][name_b]["cx"] = round(x_b, 4)
            config.protein_params[key][name_b]["cy"] = 0.0
    

    def _create_pure_membrane_system(self, system_path: str, config,
                                     molecule_imports: list) -> float:
        now          = datetime.now()
        random_int   = random.randint(100000, 999999)
        temp_dir     = os.path.join(system_path,
                                    f"temp_membrane_z-{now.strftime('%Y%m%d_%H%M%S')}-{random_int}")
        os.makedirs(temp_dir, exist_ok=True)
        try:
            membrane_params = {
                "boxx":               config.box_x,
                "boxy":               config.box_y,
                "boxz":               config.box_z,
                "box_type":           config.box_type,
                "membrane":           config.membrane_string,
                "moleculeimport":     molecule_imports,
                "solvation":          config.solvation,
                "selectedforcefield": config.selected_ff,
                "itp_input":          f"include:{self.global_info.toppar_folder_path}/{config.selected_ff}.top",
            }
            selected_molecule_builders = self._collect_selected_molecule_builders(config)
            if selected_molecule_builders:
                membrane_params["molecule_builder"] = selected_molecule_builders
            run_coby_simulation(membrane_params, "", temp_dir)
            return self._measure_membrane_z(os.path.join(temp_dir, "system.gro"), config)
        except Exception as e:
            print(f"Membrane Z measurement failed: {e}")
            return 0.0
        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def _measure_membrane_z(self, gro_path: str, config) -> float:
        try:
            u           = mda.Universe(gro_path)
            lipid_names = [entry[0] for entry in config.entries_current]
            if not lipid_names:
                return 0.0
            lipid_atoms = u.select_atoms("resname " + " ".join(lipid_names))
            if len(lipid_atoms) == 0:
                return 0.0
            z       = lipid_atoms.positions[:, 2]
            n       = max(1, int(0.1 * len(z)))
            upper_z = float(np.mean(np.sort(z)[-n:]))
            return upper_z
        except Exception as e:
            print(f"Membrane Z measurement failed: {e}")
            return 0.0

    

    def _ensure_protein_params(self, config) -> None:
        for i in range(config.n_systems):
            key = f"R{i + 1:04d}"
            if key not in config.protein_params:
                config.protein_params[key] = {}
            for pdb_name in config.protein_files:
                if pdb_name not in config.protein_params[key]:
                    config.protein_params[key][pdb_name] = {
                        "cx": 0.0, "cy": 0.0, "cz": 0.0,
                        "rx": 0.0, "ry": 0.0, "rz": 0.0,
                    }

    def _build_molecule_imports(self, config) -> list:
        imports = []
        chol = self.membrane_parser.check_chol_import(config)
        if chol == "CHOL":
            imports.append(f"file:{self.global_info.chol_file} moleculetype:{chol} params:IMPORTED")
        elif chol == "CHOL2":
            imports.append(f"file:{self.global_info.chol2_file} moleculetype:{chol} params:IMPORTED")
        if config.selected_ff == "martini_v2.2" and "WF" in config.solvation:
            imports.append(f"file:{self.global_info.wf_file} moleculetype:WF params:IMPORTED")
        return imports

    def _collect_output_pdbs(self, system_folders: list[str]) -> list[tuple[str, str]]:
        results = []
        for system_path in system_folders:
            gro_candidate = os.path.join(system_path, "eq3.gro")
            gro_path      = gro_candidate if os.path.exists(gro_candidate) \
                            else os.path.join(system_path, "system.gro")
            if not os.path.exists(gro_path):
                continue
            pdb_path = convert_gro_to_pdb(gro_path, gro_path.replace(".gro", ".pdb"))
            label    = os.path.basename(system_path)
            if pdb_path and os.path.exists(pdb_path):
                results.append((label, pdb_path))
        return results
    
    def _make_grid(self, boxx: float, boxy: float, n_copies: int) -> list[tuple[float, float]]:
        n  = math.ceil(math.sqrt(n_copies))
        dx = boxx / n
        dy = boxy / n
        return [
            (round((i + 0.5) * dx - boxx / 2, 4),
            round((j + 0.5) * dy - boxy / 2, 4))
            for i in range(n)
            for j in range(n)
        ]
    


class TopologyEditor:
    def __init__(self):
        self.global_info = GlobalInfo()
    def edit_topology(self, ff_name: str, destination: str):
        topol_file = os.path.join(destination, "topol.top")
        ff_itp = os.path.join(self.global_info.toppar_folder_path, f"{ff_name}.top")

        if not os.path.exists(ff_itp):
            print(f"{ff_itp} not found")
            return

        with open(ff_itp, 'r') as f:
            itp_lines = f.readlines()
        with open(topol_file, 'r') as f:
            topol_lines = f.readlines()[1:]

        def fix_include_path(line):
            return re.sub(
                r'(#include\s+")(.*?)(")',
                lambda m: m.group(1) + "../../toppar/" + m.group(2) + m.group(3),
                line
            )

        def collapse_blank_lines(text: str) -> str:

            return re.sub(r'\n{3,}', '\n\n', text)

        fixed_itp_lines   = [fix_include_path(line) for line in itp_lines]
        fixed_topol_lines = [fix_include_path(line) for line in topol_lines]

        merged = "".join(fixed_itp_lines) + "\n" + "".join(fixed_topol_lines)
        merged = collapse_blank_lines(merged)

        with open(topol_file, 'w') as f:
            f.write(merged)

    def replace_placeholder_title(self, topol_path: str, config) -> None:
        """
        Replaces PLACEHOLDER_TITLE in topol.top with a descriptive system title:
        {protein_filename}_{lipid1}_{upper1}_{lower1}_{lipid2}_{upper2}_{lower2}_...
        """
        # Build membrane composition string
        entries = config.entries_current if hasattr(config, "entries_current") and config.entries_current \
                else (next(iter(config.entries.values())) if isinstance(config.entries, dict)
                        else config.entries)

        composition_parts = []
        for entry in entries:
            lipid, upper, lower = entry[0], entry[1], entry[2]
            # Format ratios: drop trailing zeros for cleanliness (1.0 → 1, 0.33 → 0.33)
            upper_str = str(int(upper)) if float(upper) == int(upper) else str(upper)
            lower_str = str(int(lower)) if float(lower) == int(lower) else str(lower)
            composition_parts.append(f"{lipid}_{upper_str}_{lower_str}")

        membrane_composition = "_".join(composition_parts)

        # Build full title
        if config.protein_files:
            protein_names = "_".join(
                os.path.splitext(n)[0] for n in config.protein_files
            )
            title = f"{protein_names}_{membrane_composition}"
        else:
            title = membrane_composition

        # Replace in file
        with open(topol_path, "r") as f:
            content = f.read()

        content = content.replace("PLACEHOLDER_TITLE", title)

        with open(topol_path, "w") as f:
            f.write(content)
    def overwrite_moleculetype_line(self, itp_file: str,moleculetype_name: str = "Protein"):
        """
        Overwrites the first line after [ moleculetype ] with 'Protein'.
        
        :param itp_file: Path to the .itp file.
        """
        with open(itp_file, 'r') as f:
            lines = f.readlines()

        updated_lines = []
        inside_moleculetype = False

        for i, line in enumerate(lines):
            if re.match(r'^\s*\[\s*moleculetype\s*\]', line, re.IGNORECASE):  
                inside_moleculetype = True
                updated_lines.append(line) 
                continue
            
            if inside_moleculetype and line.strip() and not line.strip().startswith(";"):
                
                updated_lines.append(f"{moleculetype_name} 1\n")
                updated_lines.extend(lines[i + 1:])  
                break
            
            updated_lines.append(line)  

        with open(itp_file, 'w') as f:
            f.writelines(updated_lines)



    def fix_protein_includes_only(self, topol_path: str, protein_itp_names: list[str]):
        """protein_itp_names: list of stem.itp filenames e.g. ['proteinA.itp', 'proteinB.itp']"""
        with open(topol_path, 'r') as f:
            lines = f.readlines()

        updated_lines = []
        for line in lines:
            matched = False
            for itp_name in protein_itp_names:
                pattern = rf'^\s*#include\s+"[^"]*{re.escape(itp_name)}"\s*$'
                if re.match(pattern, line):
                    updated_lines.append(f'#include "../../toppar/{itp_name}"\n')
                    matched = True
                    break
            if not matched:
                updated_lines.append(line)

        with open(topol_path, 'w') as f:
            f.writelines(updated_lines)

    def fix_multiple_protein_topology(self, topol_path: str) -> None:
        with open(topol_path, "r") as f:
            lines = f.readlines()

        # ── 1. Deduplicate protein.itp includes ──────────────────────
        seen_protein_include = False
        cleaned = []
        for line in lines:
            is_protein_include = bool(re.match(r'^\s*#include\s+"[^"]*protein\.itp"', line))
            if is_protein_include:
                if not seen_protein_include:
                    cleaned.append(line)
                    seen_protein_include = True
                # else: skip duplicate
            else:
                cleaned.append(line)

        # ── 2. Merge repeated "Protein   N" lines in [ molecules ] ───
        in_molecules = False
        result = []
        protein_count = 0

        for line in cleaned:
            if re.match(r'^\s*\[\s*molecules\s*\]', line):
                in_molecules = True
                result.append(line)
                continue

            if in_molecules and re.match(r'^\s*Protein\s+\d+', line, re.IGNORECASE):
                protein_count += int(re.split(r'\s+', line.strip())[1])
                continue  # collect, don't append yet

            # flush accumulated Protein count before any non-Protein line in [ molecules ]
            if in_molecules and protein_count > 0:
                result.append(f"Protein     {protein_count}\n")
                protein_count = 0

            result.append(line)

        # flush if file ended inside [ molecules ]
        if protein_count > 0:
            result.append(f"Protein     {protein_count}\n")

        with open(topol_path, "w") as f:
            f.writelines(result)

class ForceFieldManager:
    def get_forcefield_names(self, toppar_dir: str) -> list:
        return [os.path.splitext(f)[0] for f in os.listdir(toppar_dir) if f.endswith('.top')]

    def copy_ff_folder(self, ff_name: str, destination: str):
        src_folder = os.path.join("toppar", ff_name)
        dst_folder = os.path.join(destination, ff_name)
        
        src_itp = os.path.join("toppar", f"{ff_name}.itp")
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
        src_mdp = os.path.join("mdp", ff_name, system_type)
        dst_mdp = os.path.join(destination, "mdp")
        if os.path.exists(src_mdp):
            os.makedirs(dst_mdp, exist_ok=True)
            for file in os.listdir(src_mdp):
                shutil.copy(os.path.join(src_mdp, file), os.path.join(dst_mdp, file))


def convert_gro_to_pdb(gro_file: str, pdb_file: str) -> Optional[str]:
    try:
        u = mda.Universe(gro_file)
        with mda.Writer(pdb_file) as w:
            w.write(u.atoms)
        return pdb_file
    except Exception as e:
        print(f"Gro→Pdb failed: {e}")
        return None
def run_coby_simulation(params: dict, protein_lines: Optional[str], system_path: str) -> str:
    """Combine args and run COBY"""

    box_type = params.get("box_type", params.get("boxtype", "rectangular"))

    if box_type == "rectangular":
        box = [
            float(params["boxx"]),
            float(params["boxy"]),
            float(params["boxz"]),
        ]

    elif box_type == "hexagonal":
        box = [
            float(params["boxx"]),
            float(params["boxz"]),
        ]

    elif box_type == "skewed_hexagonal":
        box = [
            float(params["boxx"]),
        ]


    else:
        raise ValueError(f"Unsupported box_type: {box_type}")


    coby_args = {
        "box_type": box_type,
        "box": box,
        "solvation": params["solvation"],
        "itp_input": params.get("itp_input"),
        "out_sys": os.path.join(system_path, "system.gro"),
        "out_top": os.path.join(system_path, "topol.top"),
    }

    if params.get("membrane"):
        coby_args["membrane"] = params["membrane"]

    if params.get("moleculeimport"):
        coby_args["molecule_import"] = params["moleculeimport"]
    if params.get("molecule_builder"):
        coby_args["molecule_builder"] = params["molecule_builder"]


    if protein_lines:
        coby_args["protein"] = protein_lines

    print(coby_args)
    COBY.COBY(**coby_args)

    return system_path
