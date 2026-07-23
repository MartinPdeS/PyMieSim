"""Source and scatterer setup section."""

from dash import dcc, html

from PyMieSim.gui.components import Card
from PyMieSim.gui.schemas import SINGLE_SCATTERER_FIELDS, SINGLE_SOURCE_FIELDS


def build_source_section():
    return _setup_card("Source Set", "single-source-type", "single-source-fields", SINGLE_SOURCE_FIELDS, "Gaussian", color="yellow", info="Configure the optical source and its wavelength, polarization, power, and aperture values.")


def build_scatterer_section():
    return _setup_card("Scatterer Set", "single-scatterer-type", "single-scatterer-fields", SINGLE_SCATTERER_FIELDS, "Sphere", color="blue", info="Configure particle geometry, material, and surrounding medium values.")


def _setup_card(title, selector_id, fields_id, choices, default, color="blue", info=""):
    return html.Details(
        id=f"{selector_id}-card",
        className=Card.classes(color=color, extra="workflow-card"),
        open=True,
        children=[
            html.Summary(children=[html.Span(title), html.Span("i", className="workflow-info-button", title=info, **{"aria-label": info})]),
            html.Div(className="workflow-card-body", children=[dcc.Dropdown(id=selector_id, className="dashboard-dropdown", options=[{"label": key, "value": key} for key in choices], value=default, clearable=False, optionHeight=38, maxHeight=200, persistence=True, persistence_type="session"), html.Div(id=fields_id, className="panel-body")]),
        ],
    )
