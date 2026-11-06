"""Page layout shared by all pages of the graphical interface.

Every page of the wizard has the same structure: a title, a sidebar with the settings and a
main area with the content. `GlobalGui.render_page` draws that frame, so the individual pages
only describe what goes into the sidebar and into the main area.
"""

import streamlit as st


class GlobalGui:
    """Building blocks that look the same on every page."""

    @staticmethod
    def render_build_button() -> bool:
        """Draw the big "BUILD!" button at the bottom of the sidebar.

        Returns:
            `True` in the moment the button is pressed, otherwise `False`.
        """
        st.sidebar.markdown("---")
        return st.sidebar.button("BUILD!", use_container_width=True, type="primary")

    @staticmethod
    def render_page(title: str, content_fn, sidebar_fn=None,
                    layout: str = "centered", sidebar_only: bool = False):
        """Draw one complete page: title, sidebar and content.

        After the page is drawn, Streamlit is told to stop, so nothing below this call is
        shown on the same page.

        Args:
            title: Page title, shown at the top and in the browser tab.
            content_fn: Function that draws the main area of the page.
            sidebar_fn: Function that draws the sidebar, or `None` for a page without one.
            layout: `"centered"` for a narrow page, `"wide"` to use the full browser width.
            sidebar_only: Draw only the title and the sidebar, without the main area.
        """
        st.markdown(
            """<style>
            div[data-testid="stBottom"] {
                position: fixed; bottom: 0; left: 0; right: 0;
                background: white; padding: 1rem 2rem;
                border-top: 1px solid #e0e0e0; z-index: 999;
            }
            </style>""",
            unsafe_allow_html=True,
        )
        st.set_page_config(page_title=title, layout=layout)
        if sidebar_fn is not None:
            with st.sidebar:
                sidebar_fn()
        st.title(title)
        if content_fn and not sidebar_only:
            st.divider()
            content_fn()
        st.stop()
