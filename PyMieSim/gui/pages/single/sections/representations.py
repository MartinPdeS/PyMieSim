"""Representation controls section."""

from dash import dcc, html

from PyMieSim.gui.components import Card
from PyMieSim.gui.defaults import DEFAULT_WORKSPACE_SETTINGS

def build_representation_section():
    return html.Section(
        id="single-representation-card",
        className=Card.classes(color="blue", extra="workflow-card representation-card"),
        children=[
            html.Div(className="representation-card-header", children=[html.Span("Representation Controls")]),
            html.Div(
                id="single-representation-body",
                className="workflow-card-body representation-card-body",
                children=[
                    html.Div(className="representation-top-row", children=[
                        html.Div(className="field-block", children=[html.Label("Representation", htmlFor="single-representation"), dcc.Dropdown(id="single-representation", className="dashboard-dropdown", options=[{"label": "S1 / S2 amplitudes", "value": "s1s2"}, {"label": "Stokes I intensity", "value": "stokes"}, {"label": "Stokes Q", "value": "stokes_q"}, {"label": "Stokes U", "value": "stokes_u"}, {"label": "Stokes V", "value": "stokes_v"}, {"label": "Scattering phase function", "value": "spf"}, {"label": "Far-field intensity", "value": "farfields"}, {"label": "Near-field E", "value": "nearfields"}, {"label": "Near-field Ex", "value": "nearfields_ex"}, {"label": "Near-field Ey", "value": "nearfields_ey"}, {"label": "Near-field Ez", "value": "nearfields_ez"}], value=DEFAULT_WORKSPACE_SETTINGS["particle_explorer"]["representation"], clearable=False, searchable=False, optionHeight=38, maxHeight=400, persistence=True, persistence_type="session")]),
                        html.Div(id="single-nearfield-mode-field", className="field-block nearfield-mode-field", children=[html.Label("Near-field values", htmlFor="single-nearfield-mode"), dcc.Dropdown(id="single-nearfield-mode", className="dashboard-dropdown", options=[{"label": "Absolute value", "value": "absolute"}, {"label": "Real part", "value": "real"}], value="absolute", clearable=False, searchable=False, optionHeight=38, maxHeight=120, persistence=True, persistence_type="session")]),
                    ]),
                    html.Div(className="field-block", children=[html.Label("Projection", htmlFor="single-projection"), dcc.Dropdown(id="single-projection", className="dashboard-dropdown", options=[{"label": "2D plot", "value": "2d"}, {"label": "Polar 1D", "value": "polar_1d"}], value=DEFAULT_WORKSPACE_SETTINGS["particle_explorer"]["projection"], clearable=False, searchable=False, optionHeight=38, maxHeight=400, persistence=True, persistence_type="session")]),
                    html.Div(className="field-block", children=[html.Label("Angular sampling", htmlFor="single-sampling"), dcc.Input(id="single-sampling", type="number", value=DEFAULT_WORKSPACE_SETTINGS["particle_explorer"]["sampling"], min=24, max=300, step=1, placeholder=str(DEFAULT_WORKSPACE_SETTINGS["particle_explorer"]["sampling"]), className="field-input", persistence=True, persistence_type="session")]),
                ],
            ),
        ],
    )
