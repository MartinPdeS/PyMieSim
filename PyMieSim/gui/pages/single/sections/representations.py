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
                    html.Div(className="field-block", children=[html.Label("Representation", htmlFor="single-representation"), dcc.Dropdown(id="single-representation", className="dashboard-dropdown", options=[{"label": "S1 / S2 amplitudes", "value": "s1s2"}, {"label": "Stokes intensity", "value": "stokes"}, {"label": "Scattering phase function", "value": "spf"}, {"label": "Far-field intensity", "value": "farfields"}], value="s1s2", clearable=False, optionHeight=34, maxHeight=200)]),
                    html.Div(className="field-block", children=[html.Label("Angular sampling", htmlFor="single-sampling"), dcc.Input(id="single-sampling", type="number", value=120, min=24, max=300, step=1, className="field-input")]),
                    html.Button("Render representation", id="run-single", n_clicks=0, className="run-button run-button-primary"),
                    html.Div(id="single-status", className="status-banner idle", children="Ready."),
                ],
            ),
        ],
    )
