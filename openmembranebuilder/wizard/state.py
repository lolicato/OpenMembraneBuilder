"""What the wizard remembers between clicks: the workflow you chose and where you are in it.

Streamlit redraws the whole page after every click. The wizard state is kept in Streamlit's
session memory, so the wizard knows which page to draw next.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
import streamlit as st


class WizardStep(Enum):
    """The pages of the wizard.

    Attributes:
        SELECT_PATHWAY: The welcome page, where you choose what to build.
        MARTINIZE: Preparing the proteins (upload coarse-grained files or martinize an
            atomistic structure).
        MEMBRANE: Composition, box, solvation, protein placement and preview.
        RESULTS: Download of the zip archive and 3D view of the built systems.
    """

    SELECT_PATHWAY = auto()
    MARTINIZE      = auto()
    MEMBRANE       = auto()
    RESULTS        = auto()


@dataclass(frozen=True)
class Pathway:
    """One workflow offered on the welcome page.

    Attributes:
        label: Text shown on the welcome page.
        steps: The pages this workflow goes through, in order.
        has_protein: The workflow includes proteins.
        protein_only: Proteins in water, without a membrane.
    """

    label:        str
    steps:        list[WizardStep]
    has_protein:  bool = False
    protein_only: bool = False


PATHWAYS: list[Pathway] = [
    Pathway(
        label="Membrane only (no protein)",
        steps=[WizardStep.MEMBRANE],
        has_protein=False,
    ),
    Pathway(
        label="Protein + Membrane",
        steps=[WizardStep.MARTINIZE, WizardStep.MEMBRANE],
        has_protein=True,
    ),
    Pathway(
        label="Protein only (no membrane)",
        steps=[WizardStep.MARTINIZE, WizardStep.MEMBRANE],
        has_protein=True,
        protein_only=True,
    ),
]


@dataclass
class WizardState:
    """Where the user is in the wizard.

    Attributes:
        step: The page shown at the moment.
        pathway: The chosen workflow, `None` before the welcome page is left.
        step_index: Position of the current page within the workflow.
        temp_folder: Folder for the uploads of this session.
        completed: The build has finished.
        locked: A workflow has been started; the sidebar then offers "Change workflow".
    """

    step:                 WizardStep   = WizardStep.SELECT_PATHWAY
    pathway:              Pathway|None = None
    step_index:           int          = 0
    temp_folder:          str          = field(default_factory=lambda: "")
    completed:            bool         = False
    locked:               bool         = False

    @property
    def has_protein(self) -> bool:
        """`True` when the chosen workflow contains proteins."""
        return self.pathway.has_protein if self.pathway else False

    @property
    def needs_martinize(self) -> bool:
        """`True` when the chosen workflow includes the protein-preparation step."""
        return self.pathway is not None and WizardStep.MARTINIZE in self.pathway.steps


def get_state() -> WizardState:
    """The wizard state of this browser session, created on first use.

    Returns:
        The `WizardState` kept in Streamlit's session memory.
    """
    if "wizard_state" not in st.session_state:
        st.session_state.wizard_state = WizardState()
    return st.session_state.wizard_state