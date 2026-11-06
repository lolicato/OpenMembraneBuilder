"""Allows `python -m openmembranebuilder gui|cli ...`, which does exactly the same as `omb gui|cli ...`.

Useful when the `omb` command is not on the PATH (for example in some cluster environments).
"""
import sys

from openmembranebuilder.cli import main

sys.exit(main())
