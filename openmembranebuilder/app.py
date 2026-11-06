"""The page Streamlit opens when the graphical interface starts.

`omb gui` asks Streamlit to run this file. All it does is hand over to the wizard
(`openmembranebuilder.wizard.page.render_wizard`), which decides which page to show.
"""
from openmembranebuilder.wizard.page import render_wizard

render_wizard()
