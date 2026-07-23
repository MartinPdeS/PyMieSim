"""Representation controls section."""

from dash import dcc, html

from PyMieSim.gui.components import Card

def build_representation_section():
    return html.Section(
        id="single-representation-card",
        className=Card.classes(color="blue", extra="workflow-card representation-card"),
        children=[
            html.Div(className="representation-card-header", children=[html.Span("Representation Controls")]),
            html.Div(
                className="workflow-card-body representation-card-body",
                children=[
                    html.Div(className="field-block", children=[html.Label("Representation", htmlFor="single-representation"), dcc.Dropdown(id="single-representation", className="dashboard-dropdown", options=[{"label": "S1 / S2 amplitudes", "value": "s1s2"}, {"label": "Stokes I intensity", "value": "stokes"}, {"label": "Stokes Q", "value": "stokes_q"}, {"label": "Stokes U", "value": "stokes_u"}, {"label": "Stokes V", "value": "stokes_v"}, {"label": "Scattering phase function", "value": "spf"}, {"label": "Far-field intensity", "value": "farfields"}], value="s1s2", clearable=False, optionHeight=38, maxHeight=400, persistence=True, persistence_type="session")]),
                    html.Div(className="field-block", children=[html.Label("Projection", htmlFor="single-projection"), dcc.Dropdown(id="single-projection", className="dashboard-dropdown", options=[{"label": "2D heatmap", "value": "2d"}, {"label": "3D surface", "value": "3d"}], value="2d", clearable=False, optionHeight=38, maxHeight=400, persistence=True, persistence_type="session")]),
                    html.Div(className="field-block", children=[html.Label("Angular sampling", htmlFor="single-sampling"), dcc.Input(id="single-sampling", type="number", value=120, min=24, max=300, step=1, placeholder="120", className="field-input", persistence=True, persistence_type="session")]),
                ],
            ),
        ],
    )
