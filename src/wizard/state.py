from dataclasses import dataclass, field
from enum import Enum, auto
import streamlit as st


class WizardStep(Enum):
    SELECT_PATHWAY = auto()
    MARTINIZE      = auto()
    MEMBRANE       = auto()
    RESULTS        = auto()


@dataclass(frozen=True)
class Pathway:
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
    step:                 WizardStep   = WizardStep.SELECT_PATHWAY
    pathway:              Pathway|None = None
    step_index:           int          = 0
    temp_folder:          str          = field(default_factory=lambda: "")
    completed:            bool         = False
    locked:               bool         = False

    @property
    def has_protein(self) -> bool:
        return self.pathway.has_protein if self.pathway else False

    @property
    def needs_martinize(self) -> bool:
        return self.pathway is not None and WizardStep.MARTINIZE in self.pathway.steps


def get_state() -> WizardState:
    if "wizard_state" not in st.session_state:
        st.session_state.wizard_state = WizardState()
    return st.session_state.wizard_state