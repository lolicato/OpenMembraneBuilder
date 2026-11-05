"""Kept for the old commands; prefer `omb gui` and `omb cli config.json`.

    streamlit run app.py                 (= omb gui)
    python app.py --no-gui config.json   (= omb cli config.json)
"""
import sys


def main():
    args = sys.argv[1:]
    if "--no-gui" in args:
        idx = args.index("--no-gui")
        config_path = args[idx + 1] if idx + 1 < len(args) else None
        if not config_path:
            print("Error: --no-gui requires a config path")
            sys.exit(1)
        from openmembranebuilder.wizard.runner import run_from_config
        run_from_config(config_path)
    else:
        from openmembranebuilder.wizard.page import render_wizard
        render_wizard()


main()
