"""Parameter Sweep page composition."""

from __future__ import annotations

from dash import dcc, html

from PyMieSim.gui.components import HeaderCard
from PyMieSim.gui.layout import PLOT_CONFIG, _plot_options_card, _x_axis_card

from .sections import build_detector_section, build_scatterer_section, build_source_section


def build_experiment_page(default_measure_options: list[str], plot_settings: dict | None = None):
    """Build the isolated Parameter Sweep workspace."""
    settings = (plot_settings or {}).get("parameter_sweep", {})
    return html.Div(
        className="tab-content-stack experiment-tab-content",
        children=[
            HeaderCard(
                "Parameter Sweep Lab",
                "Configure source, scatterer, and detector sets, then run parameter sweeps through the compiled engine directly from Dash.",
                [
                    ("01", "Configure source", "Choose the source family and sweep its optical parameters.", "yellow"),
                    ("02", "Configure scatterer", "Set particle geometry, material, and medium values.", "blue"),
                    ("03", "Configure detector", "Select collection geometry and compute the sweep response.", "orange"),
                ],
                color="green",
            ).render(),
            html.Section(
                id="configure",
                className="workspace",
                children=[
                    html.Section(
                        className="control-column",
                        children=[html.Div(className="set-panel-grid", children=[build_source_section(), build_scatterer_section(), build_detector_section()])],
                    ),
                    html.Section(
                        className="result-column",
                        children=[
                            html.Section(className="panel graph-panel", children=[dcc.Loading(id="result-graph-loading", type="circle", color="#4f8df7", custom_spinner=html.Div("Computing…", className="plot-computing-indicator"), delay_show=150, delay_hide=150, children=html.Div(className="plot-loading-target", children=[dcc.Graph(id="result-graph", config=PLOT_CONFIG), html.Div(id="experiment-computation-status", style={"display": "none"})]))]),
                            _x_axis_card(default_measure_options),
                            html.Div(id="experiment-plot-options-container", children=[_plot_options_card("experiment", settings)]),
                            html.Div(className="export-actions", children=[html.Button("Export CSV", id="export-csv", n_clicks=0, className="run-button export-button")]),
                        ],
                    ),
                ],
            ),
        ],
    )
