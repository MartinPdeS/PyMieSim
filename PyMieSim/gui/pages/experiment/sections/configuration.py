"""Source, scatterer, and detector configuration cards."""

from dash import dcc, html

from PyMieSim.gui.components import Card
from PyMieSim.gui.schemas import SECTION_FIELDS


def build_source_section():
    return _build_section("Source Set", "source-type", "source-fields", SECTION_FIELDS["source"], "PlaneWaveSet", color="yellow", info="Configure the optical source and its wavelength, polarization, and amplitude values.")


def build_scatterer_section():
    return _build_section("Scatterer Set", "scatterer-type", "scatterer-fields", SECTION_FIELDS["scatterer"], "SphereSet", color="blue", info="Configure particle geometry, material, and surrounding medium values.")


def build_detector_section():
    return _build_section("Detector Set", "detector-type", "detector-fields", SECTION_FIELDS["detector"], "PhotodiodeSet", detector=True, color="cyan", info="Configure detector collection geometry, sampling, offsets, and mode settings.")


def _build_section(title, selector_id, fields_id, choices, default, detector=False, color="blue", info=""):
    options = [{"label": "No detector" if detector and key == "None" else key.replace("Set", ""), "value": key} for key in choices]
    return html.Details(
        id=f"{selector_id}-card",
        className=Card.classes(color=color, extra="workflow-card"),
        open=False,
        children=[
            html.Summary(children=[html.Span(title), html.Span("i", className="workflow-info-button", title=info, **{"aria-label": info})]),
            html.Div(className="workflow-card-body", children=[dcc.Dropdown(id=selector_id, className="dashboard-dropdown", options=options, value=default, clearable=False, optionHeight=38, maxHeight=200, persistence="parameter-sweep-defaults-v2", persistence_type="session"), html.Div(id=fields_id, className="panel-body")]),
        ],
    )
