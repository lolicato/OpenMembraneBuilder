# 🧾 OpenMembraneBuilder-v1.2.0-dev

## Added

* Installable with `pip` (`pip install openmembranebuilder-1.2.0.dev0.tar.gz`, or `pip install .`)
* The `omb` command:
    ```bash
    omb gui                 # graphical app (extra options are passed to Streamlit)
    omb cli config.json     # build from a configuration
    omb --version
    ```
* Protein insertion page: system view (x–z, y–z, x–y) and box size
* Protein names (`<name>.pdb`, `<name>.itp`, moleculetype `<name>`), copies moved together or independently

## Changed

* The code lives in the `openmembranebuilder` package (formerly `src/`), with the
  force-field resources inside it, so they are installed with the package
* Outputs, downloads and temporary uploads go to the folder `omb` is run from (or `$OMB_WORKDIR`)
* `streamlit run app.py` and `python app.py --no-gui config.json` still work

---

# 🧾 OpenMembraneBuilder-v1.0.3-dev

## Added

* Command-line interface (CLI) support
    * Run builds without GUI using:
      ```bash
      python app.py --no-gui config.json
      ```
    * JSON configuration file is generated automatically for reuse and scripting

* JSON-based workflow
    * Full build configuration is exported to `config.json`
    * JSON is stored inside the output project folder
    * Included in the final zipped archive for reproducibility

* Automatic cleanup of temporary upload folders

    * `temp_uploads-*` directories are now removed after each build
    * Cleanup is handled via a `try / except / finally` block to ensure reliability

* Support for cardiolipin lipids (TMCL, TOCL)

---

## Changed

* Output structure
    * All builds are now created inside an `outputs/` directory
    * Each build is self-contained:
      ```
      outputs/OpenMembraneBuilder-xxxx/
      ```

* Input file handling

    * Uploaded PDB and ITP files are now copied into the output project folder (`user_inputs/`)
    * Ensures all required inputs are included in the downloadable ZIP

* ZIP packaging behavior
    * ZIP now includes:
        * simulation data
        * topology files
        * `config.json`
    * Packaging happens after JSON creation

* Protein input handling
    * Uploaded PDB and ITP files are copied into the project folder
    * Paths are updated to be portable and reusable in CLI mode

* Solvation handling
    * Salt molarity is now consistently managed via `Config`
    * Eliminated inconsistencies between string and numeric representations

---

## Fixed

* Removed dependency on `st.session_state` outside GUI
    * All core logic now relies on `Config` object
    * Enables proper CLI execution


---

## Internal

* Refactored builders and inserter modules
    * Fully decoupled from Streamlit
    * Now operate purely on `Config`
