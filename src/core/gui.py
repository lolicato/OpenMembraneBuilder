import streamlit as st


class GlobalGui:

    @staticmethod
    def render_build_button() -> bool:
        st.sidebar.markdown("---")
        return st.sidebar.button("BUILD!", use_container_width=True, type="primary")

    @staticmethod
    def render_page(title: str, content_fn, sidebar_fn=None,
                    layout: str = "centered", sidebar_only: bool = False):
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
