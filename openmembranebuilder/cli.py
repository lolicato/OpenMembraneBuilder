"""The `omb` command, the single entry point of OpenMembraneBuilder.

```text
omb gui [--port 8501] [streamlit options]   start the graphical interface in the browser
omb gui --no-run                            the same, without the "Run simulations" button
omb cli config.json                         rebuild the systems of a saved config.json (no GUI)
omb cli template.json --pdb 1UBQ            the same, with another protein (prepared with martinize2)
omb --version                               print the installed version
```

`omb gui` starts Streamlit with the wizard page; any option that `omb` itself does not know
(for example `--server.headless true`) is passed on to Streamlit unchanged.

`omb gui --no-run` hides the button that starts the simulations after a build – use it when
the app is offered to others, e.g. on a web server, so that nobody can run programs on it.

`omb cli` reads a `config.json` written by an earlier build and builds the same systems again,
without opening a browser – handy for scripts, clusters and batch work.

Everything the program writes (the `output/` and `downloads/` folders) goes to the folder the
command is started from, or to the folder named by the `OMB_WORKDIR` environment variable.
"""
import argparse
import os
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version


def _version() -> str:
    """The installed version of OpenMembraneBuilder, as recorded by pip.

    Returns:
        The version text, e.g. `"2.0.0"`, or a short note when the package is used
        without being installed (for example straight from a source folder).
    """
    try:
        return version("OpenMembraneBuilder")
    except PackageNotFoundError:
        return "unknown (not installed)"


def main(argv=None) -> int:
    """Run the `omb` command.

    Reads the command line, then either starts the graphical interface (`gui`) or rebuilds
    the systems of a configuration file (`cli`).

    Args:
        argv: The command-line words after `omb`. Leave empty to use the real command line;
            tests pass a list such as `["cli", "config.json"]`.

    Returns:
        The exit code: `0` when everything went well, `1` when a protein could not be
        prepared, `3` when the configuration file does
        not exist, otherwise the exit code of Streamlit.
    """
    parser = argparse.ArgumentParser(prog="omb", description="OpenMembraneBuilder: Martini membrane and membrane–protein systems with COBY")
    parser.add_argument("--version", action="version", version=f"OpenMembraneBuilder {_version()}")
    sub = parser.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gui", help="start the graphical interface (Streamlit) in the browser")
    g.add_argument("--port", type=int, default=None, help="port of the web server (default 8501)")
    g.add_argument("--no-run", action="store_true",
                   help="hide the 'Run simulations' button (e.g. when the app is offered on a web server)")
    c = sub.add_parser("cli", help="build the systems of a config.json without the GUI")
    c.add_argument("config", help="config.json written by an earlier build (or a template)")
    c.add_argument("--pdb", action="append", metavar="ID_OR_FILE",
                   help="use this structure for the protein: a PDB code (downloaded from RCSB) or an "
                        "atomistic .pdb file, prepared with the martinize2 settings of the config; "
                        "repeat for several proteins")
    c.add_argument("--martinize", action="store_true",
                   help="prepare the proteins again with martinize2, even if finished files are recorded")
    args, extra = parser.parse_known_args(argv)

    if args.cmd == "gui":
        app = os.path.join(os.path.dirname(os.path.abspath(__file__)), "app.py")
        cmd = [sys.executable, "-m", "streamlit", "run", app]
        if args.port:
            cmd += ["--server.port", str(args.port)]
        env = dict(os.environ)
        if args.no_run:  # read by the results page (wizard/run_simulations.py)
            env["OMB_NO_RUN"] = "1"
        return subprocess.call(cmd + extra, env=env)  # extra: any other `streamlit run` option

    if extra:
        parser.error(f"unrecognized arguments: {' '.join(extra)}")
    if not os.path.exists(args.config):
        print(f"Error: {args.config} not found", file=sys.stderr)
        return 3
    from openmembranebuilder.wizard.runner import run_from_config
    try:
        run_from_config(args.config, pdbs=args.pdb, martinize=args.martinize)
    except (RuntimeError, ValueError) as e:      # e.g. a protein that could not be prepared
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
