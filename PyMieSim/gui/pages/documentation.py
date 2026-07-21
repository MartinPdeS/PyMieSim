"""Documentation hub page composition."""

from dash import html

from PyMieSim.gui.components import Card


def build_documentation_page():
    """Build the documentation index page."""
    topics = (
        ("Experiment workflow", "Use source sets, scatterer sets, and detector sets to build sweeps. The X-axis selector is inferred from fields containing multiple values."),
        ("Single workflow", "Use one source and one scatterer to inspect angular amplitudes, polarization, phase functions, or far-field intensity."),
        ("Field syntax", "Quantity inputs accept scalar values, comma-separated lists, and start:end:count expressions such as 400:1400:8."),
        ("Supported models", "Gaussian and plane-wave sources are available alongside spheres, infinite cylinders, and core-shell scatterers."),
    )
    return html.Div(
        className="page-content-stack",
        children=[
            html.Section(className="page-hero documentation-page-hero", children=[html.P("PyMieSim reference", className="eyebrow"), html.H1("Documentation"), html.P("A practical map of the objects, workflows, and input conventions used throughout the dashboard.", className="hero-text")]),
            html.Div(className="documentation-grid", children=[_topic_card(title, description) for title, description in topics]),
            html.Section(className=Card.classes(color="blue", extra="panel documentation-note"), children=[html.Div(className="card-header panel-header", children=[html.H2("Where to start")]), html.Div(className="card-body", children=[html.P("Choose Experiment for parameter studies and detector coupling. Choose Single for representation plots and physical intuition."), html.A("Open the Experiment workspace →", href="/experiment", className="inline-action")])]),
        ],
    )


def _topic_card(title: str, description: str):
    return html.Section(className=Card.classes(color="blue", extra="panel documentation-card"), children=[html.Div(className="card-header panel-header", children=[html.H2(title)]), html.Div(className="card-body", children=[html.P(description)])])
