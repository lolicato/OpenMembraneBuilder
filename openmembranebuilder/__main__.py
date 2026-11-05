"""`python -m openmembranebuilder gui|cli ...` is the same as `omb gui|cli ...`."""
import sys

from openmembranebuilder.cli import main

sys.exit(main())
