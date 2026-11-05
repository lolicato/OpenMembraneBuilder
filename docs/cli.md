# CLI Usage

OpenMembraneBuilder can run without the GUI using a saved JSON configuration.

```bash
omb cli config.json
```

---

## 🎯 Why use the CLI

The CLI mode is useful for:

- Repeating a previous build  
- Quickly swapping proteins  
- Automating workflows  
- Running batch simulations  

---

## 🔁 Typical Workflow

### 1. Create a configuration via GUI

```bash
omb gui
```

Build your system once. This generates a project folder with a `config.json`.

---

### 2. Re-run using CLI

```bash
omb cli outputs/OpenMembraneBuilder-xxxx/config.json
```

---

## ✏️ Editing configurations

You can modify `config.json` to:

- Change protein files  
- Adjust membrane composition  
- Modify box size or number of systems  

This enables fast iteration without reopening the GUI.

---

## 📦 Output

The CLI produces the same outputs as the GUI:

- Project folder in `outputs/`  
- ZIP archive in `downloads/`  
- Full reproducibility via `config.json`  
