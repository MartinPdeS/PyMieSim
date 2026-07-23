"""Single-particle page composition."""

from dash import dcc, html

from PyMieSim.gui.components import HeaderCard
from PyMieSim.gui.layout import PLOT_CONFIG, _plot_options_card

from .sections import build_representation_section, build_scatterer_section, build_source_section


def build_single_page(plot_settings: dict | None = None):
    """Build the single page from independent setup and plot sections."""
    return html.Div(
        className="tab-content-stack single-tab-content",
        children=[
            HeaderCard(
                "Particle Explorer",
                "Inspect angular scattering, polarization, and field patterns for one source–scatterer setup.",
                [
                    ("01", "Configure source", "Choose a source and set its optical parameters.", "yellow"),
                    ("02", "Configure scatterer", "Define the particle geometry and refractive indices.", "blue"),
                    ("03", "Choose representation", "Inspect amplitudes, polarization, phase functions, or far fields.", "orange"),
                ],
                color="green",
            ).render(),
            html.Section(className="single-workspace", children=[html.Section(className="control-column", children=[build_source_section(), build_scatterer_section()]), html.Section(className="result-column", children=[html.Section(className="panel graph-panel single-graph-panel", children=[dcc.Loading(id="single-graph-loading", type="circle", color="#4f8df7", custom_spinner=html.Div("Computing…", className="plot-computing-indicator"), delay_show=150, delay_hide=150, children=html.Div(className="plot-loading-target", children=[dcc.Graph(id="single-graph", config=PLOT_CONFIG), html.Div(id="single-computation-status", style={"display": "none"})]))]), build_representation_section(), html.Div(id="single-plot-options-container", children=[_plot_options_card("single", plot_settings or {})]), html.Div(className="export-actions", children=[html.Button("Export CSV", id="export-single-csv", n_clicks=0, className="run-button export-button")])])]),
        ],
    )
