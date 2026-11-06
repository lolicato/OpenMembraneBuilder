"""The wizard: the step-by-step workflow of the graphical interface.

- `state` – what the wizard remembers (chosen workflow, current page).
- `navigation` – moving forward and back.
- `page` – the welcome page, the hand-over between builders and saving `config.json`.
- `gui` – the results page (download and 3D view).
- `run_simulations` – the "Run simulations" button of the results page.
- `runner` – the same build without the interface, from a `config.json` (`omb cli`).
"""

