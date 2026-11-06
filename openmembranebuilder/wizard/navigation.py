"""Moving forward and backward through the wizard.

A workflow ("pathway") is a short list of steps, for example *protein preparation → membrane*.
These two functions move the wizard one step along that list, or back, and bring it to the
results page after the last step.
"""

from openmembranebuilder.wizard.state import WizardState, WizardStep


def go_forward(state: WizardState) -> None:
    """Move the wizard to the next step of the chosen workflow.

    From the welcome page this is the first step of the workflow; after the last step it is
    the results page.

    Args:
        state: The wizard's state (see `WizardState`), changed in place.
    """
    if state.step == WizardStep.SELECT_PATHWAY:
        state.step_index = 0
        _sync_step(state)
        return

    state.step_index += 1
    if state.step_index >= len(state.pathway.steps):
        state.step = WizardStep.RESULTS
    else:
        _sync_step(state)


def go_back(state: WizardState) -> None:
    """Move the wizard one step back.

    From the first step of a workflow this returns to the welcome page, where another
    workflow can be chosen.

    Args:
        state: The wizard's state (see `WizardState`), changed in place.
    """
    if state.step_index > 0:
        state.step_index -= 1
        _sync_step(state)
        return

    # Already on the first step of the pathway -> return to the pathway selection
    # (previously this was a silent no-op and the user was stuck).
    state.step = WizardStep.SELECT_PATHWAY
    state.step_index = 0
    state.locked = False


def _sync_step(state: WizardState) -> None:
    """Set the current page to the step at the current position of the workflow.

    Args:
        state: The wizard's state, changed in place.
    """
    state.step = state.pathway.steps[state.step_index]
