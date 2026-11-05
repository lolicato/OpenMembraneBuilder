# GUI Usage

OpenMembraneBuilder's graphical app is started with:

```bash
omb gui
```

---

## 🧭 What you can configure

The GUI allows you to configure:

- System type (membrane or membrane + protein)  
- Box size  
- Martini force field  
- Membrane composition  
- Solvation and salt molarity  
- Number of replicas (systems)  
- Protein input files  
- Protein placement and rotation  

---

## 🚀 Build Process

Press the **BUILD** button to generate your systems.

This will create a project folder:

```
outputs/<project_name>/
```

and automatically generate a downloadable ZIP archive.

---

## 🧬 Protein Input Files

For protein-containing systems, upload:

- A coarse-grained `.pdb` file  
- The corresponding `.itp` topology file  

These files are copied into:

```
user_inputs/
```

inside the output project folder, ensuring the build is fully reproducible.

---

## 📦 Output Summary

After building, you will have:

- A complete project folder inside `outputs/`  
- A ZIP archive in `downloads/`  
- All inputs and configuration stored for reuse  
