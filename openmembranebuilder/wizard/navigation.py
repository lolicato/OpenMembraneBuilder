from openmembranebuilder.wizard.state import WizardState, WizardStep


def go_forward(state: WizardState) -> None:
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
    state.step = state.pathway.steps[state.step_index]
