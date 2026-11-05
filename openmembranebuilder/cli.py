"""The `omb` command.

    omb gui [--port 8501] [streamlit options]   start the graphical interface in the browser
    omb cli config.json                         build the systems of a saved config.json (no GUI)
    omb --version
"""
import argparse
import os
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version


def _version() -> str:
    try:
        return version("OpenMembraneBuilder")
    except PackageNotFoundError:
        return "unknown (not installed)"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="omb", description="OpenMembraneBuilder: Martini membrane and membrane–protein systems with COBY")
    parser.add_argument("--version", action="version", version=f"OpenMembraneBuilder {_version()}")
    sub = parser.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gui", help="start the graphical interface (Streamlit) in the browser")
    g.add_argument("--port", type=int, default=None, help="port of the web server (default 8501)")
    c = sub.add_parser("cli", help="build the systems of a config.json without the GUI")
    c.add_argument("config", help="config.json written by an earlier build")
    args, extra = parser.parse_known_args(argv)

    if args.cmd == "gui":
        app = os.path.join(os.path.dirname(os.path.abspath(__file__)), "app.py")
        cmd = [sys.executable, "-m", "streamlit", "run", app]
        if args.port:
            cmd += ["--server.port", str(args.port)]
        return subprocess.call(cmd + extra)  # extra: any other `streamlit run` option

    if extra:
        parser.error(f"unrecognized arguments: {' '.join(extra)}")
    if not os.path.exists(args.config):
        print(f"Error: {args.config} not found", file=sys.stderr)
        return 3
    from openmembranebuilder.wizard.runner import run_from_config
    run_from_config(args.config)
    return 0


if __name__ == "__main__":
    sys.exit(main())
