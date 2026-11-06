"""Starting the simulations from the results page (the "Run" button).

After a build, the results page offers **▶ Run simulations**. It starts `systems/run.sh` of
the build – energy minimisation, index groups, equilibration and production for every system –
on the computer where the app runs. The script runs in the background: you can close the page,
the simulations continue. Its output goes to `systems/run.log`, its process number to
`systems/run.pid`, so the page shows the progress and can stop it again.

The button is not shown when the app is started with `omb gui --no-run` (or with the
environment variable `OMB_NO_RUN=1`), for example on a web server where visitors must not run
anything on the server.
"""

from __future__ import annotations

import os
import shutil
import signal
import subprocess

import streamlit as st

# Processes started by this app, by process number. Keeping them lets the app collect a
# finished process, so it is not mistaken for a running one.
_STARTED: dict[int, subprocess.Popen] = {}

#: Programs `run.sh` needs.
REQUIRED_PROGRAMS = ("gmx", "jq")


def run_allowed() -> bool:
    """Whether the Run button may be offered.

    It is switched off by `omb gui --no-run`, which sets the environment variable
    `OMB_NO_RUN=1` for the app. The variable can also be set by hand, e.g. when the app is
    started with `streamlit run app.py` on a server.

    Returns:
        `False` when `OMB_NO_RUN` is `1`, `true`, `yes` or `on`, otherwise `True`.
    """
    return os.environ.get("OMB_NO_RUN", "").strip().lower() not in ("1", "true", "yes", "on")


def missing_programs() -> list[str]:
    """The programs `run.sh` needs that cannot be found on this computer.

    Returns:
        The missing names among `gmx` and `jq` (empty when both are found).
    """
    return [p for p in REQUIRED_PROGRAMS if shutil.which(p) is None]


def _paths(output_folder: str) -> tuple[str, str, str, str]:
    """The `systems` folder of a build, and its `run.sh`, `run.log` and `run.pid`."""
    systems = os.path.join(output_folder, "systems")
    return (systems, os.path.join(systems, "run.sh"),
            os.path.join(systems, "run.log"), os.path.join(systems, "run.pid"))


def running_pid(pid_file: str) -> int | None:
    """The process number of the simulations, if they are still running.

    Args:
        pid_file: The `run.pid` written when the simulations were started.

    Returns:
        The process number, or `None` when nothing was started or it has finished.
    """
    try:
        pid = int(open(pid_file).read().strip())
    except (OSError, ValueError):
        return None
    proc = _STARTED.get(pid)
    if proc is not None:  # started by this app: ask it directly
        return pid if proc.poll() is None else None
    try:  # started earlier (e.g. before the app was restarted)
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return None
    return pid


def start_run(output_folder: str) -> int:
    """Start `run.sh` in the `systems` folder of a build, in the background.

    The script keeps running when the page is closed. Its output is added to
    `systems/run.log`; its process number is written to `systems/run.pid`.

    Args:
        output_folder: The output folder of the build.

    Returns:
        The process number of the started script.
    """
    systems, script, log_file, pid_file = _paths(output_folder)
    with open(log_file, "ab") as log:
        proc = subprocess.Popen(
            ["bash", os.path.basename(script)],
            cwd=systems,                       # run.sh runs every system folder in here
            stdout=log,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=True,            # own process group: survives the page, can be stopped as a whole
        )
    _STARTED[proc.pid] = proc
    with open(pid_file, "w") as f:
        f.write(str(proc.pid))
    return proc.pid


def stop_run(pid: int) -> None:
    """Stop the simulations: `run.sh` and the GROMACS programs it started.

    Args:
        pid: The process number of `run.sh`.
    """
    try:
        os.killpg(pid, signal.SIGTERM)       # the whole group: run.sh and gmx
    except (AttributeError, OSError):         # no process groups (Windows) or already gone
        proc = _STARTED.get(pid)
        if proc is not None:
            proc.terminate()


def _tail(path: str, n: int = 40) -> str:
    """The last `n` lines of a text file (empty when it does not exist)."""
    try:
        with open(path, errors="replace") as f:
            return "".join(f.readlines()[-n:])
    except OSError:
        return ""


def render_run_controls(output_folder: str | None) -> None:
    """The "Run simulations" part of the results page.

    Shows **▶ Run simulations** while nothing runs, and the state, the latest lines of
    `run.log` and **⏹ Stop** while the simulations run. Nothing is shown when running is
    switched off (`omb gui --no-run`), or when the build has no `systems/run.sh`.

    The button is disabled, with an explanation, when `gmx` or `jq` cannot be found.

    Args:
        output_folder: The output folder of the build that was just made.
    """
    if not run_allowed() or not output_folder:
        return
    systems, script, log_file, pid_file = _paths(output_folder)
    if not os.path.exists(script):
        return

    st.subheader("Run the simulations")
    st.caption(
        f"Starts `run.sh` in `{systems}`: energy minimisation, equilibration (eq1–eq3) and "
        "production for every system, on this computer. It keeps running when you close the page."
    )

    pid = running_pid(pid_file)
    if pid is not None:
        st.info(f"⏳ The simulations are running (process {pid}).")
        left, right, _ = st.columns([1, 1, 4])
        if left.button("🔄 Refresh", use_container_width=True, key="run_refresh"):
            st.rerun()
        if right.button("⏹ Stop", use_container_width=True, key="run_stop"):
            stop_run(pid)
            st.rerun()
    else:
        missing = missing_programs()
        if missing:
            st.warning(
                f"Cannot run here: {', '.join(missing)} not found. Make GROMACS (`gmx`) and "
                "`jq` available, then restart `omb gui` – or run `systems/run.sh` yourself."
            )
        elif os.path.exists(log_file):
            st.success("The last run has ended – see the log below.")
        if st.button("▶ Run simulations", type="primary", disabled=bool(missing), key="run_start"):
            start_run(output_folder)
            st.rerun()

    log = _tail(log_file)
    if log:
        with st.expander("run.log (last lines)", expanded=pid is not None):
            st.code(log, language="text")
