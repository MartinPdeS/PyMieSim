"""Home page composition."""

from dash import html

from PyMieSim.gui.components import Card

def build_home_page():
    """Build the landing page and its workflow cards."""
    return html.Div(
        className="page-content-stack",
        children=[
            html.Section(
                className="page-hero",
                children=[
                    html.P("PyMieSim", className="eyebrow"),
                    html.H1("Optical scattering, clearly organized."),
                    html.P(
                        "Explore parameter sweeps in the Experiment workspace or inspect a single particle through its physical representations.",
                        className="hero-text",
                    ),
                ],
            ),
            _citation_card(),
            html.Div(
                className="home-capability-grid",
                children=[
                    _capability_card(
                        "Single module",
                        "Inspect one optical setup through angular scattering, polarization, phase functions, and field representations.",
                        ["Configure source", "Configure scatterer", "Render representations"],
                        "Open Single module",
                        "/single",
                    ),
                    _capability_card(
                        "Experiment",
                        "Run source, scatterer, and detector sets across parameter sweeps and export structured results for analysis.",
                        ["Configure source", "Configure scatterer and detector", "Run and export results"],
                        "Open experiment",
                        "/experiment",
                    ),
                    _capability_card(
                        "Documentation",
                        "Learn the PyMieSimX vocabulary, field syntax, supported objects, and recommended workflows.",
                        ["Learn the model vocabulary", "Explore field syntax", "Follow recommended workflows"],
                        "Read documentation",
                        "/documentation",
                    ),
                ],
            ),
            _metrics_card(),
        ],
    )


def _citation_card():
    """Build the support and citation panel shown on the landing page."""
    return html.Section(
        className=Card.classes(color="blue", extra="home-info-card"),
        children=[
            html.Div("Support, citation, and lab", className="home-section-header"),
            html.Div(
                className="home-info-body",
                children=[
                    html.P(
                        "Support PyMieSim development, cite the underlying work in publications, and learn more about the project.",
                        className="home-info-copy",
                    ),
                    html.Div(
                        className="home-button-row",
                        children=[
                            html.A("Support Developer", href="https://github.com/MartinPdeS/PyMieSim", target="_blank", rel="noopener noreferrer", className="home-button home-button-primary"),
                            html.A("Citing this work", href="https://opg.optica.org/optcon/fulltext.cfm?uri=optcon-2-3-520", target="_blank", rel="noopener noreferrer", className="home-button home-button-outline"),
                            html.A("PyMieSim documentation", href="https://martinpdes.github.io/PyMieSim/", target="_blank", rel="noopener noreferrer", className="home-button home-button-info"),
                            html.A("GitHub repository", href="https://github.com/MartinPdeS/PyMieSim", target="_blank", rel="noopener noreferrer", className="home-button home-button-muted"),
                        ],
                    ),
                ],
            ),
        ],
    )


def _capability_card(title: str, description: str, steps: list[str], button_text: str, href: str):
    return html.Section(
        className=Card.classes(color="blue", extra="home-capability-card"),
        children=[
            html.Div(title, className="home-section-header"),
            html.Div(
                className="home-capability-body",
                children=[
                    html.P(description),
                    html.Div(
                        [html.Div([html.Span(str(index), className="home-step-number"), html.Span(step)]) for index, step in enumerate(steps, start=1)],
                        className="home-step-list",
                    ),
                    html.A(button_text, href=href, className="home-workflow-button"),
                ],
            ),
        ],
    )


def _metrics_card():
    return html.Section(
        className=Card.classes(color="blue", extra="home-info-card home-metrics-card"),
        children=[
            html.Div("PyMieSim at a glance", className="home-section-header"),
            html.Div(
                className="home-metrics-grid",
                children=[
                    html.Div([html.Strong("2", className="home-metric-value"), html.Span("analysis workspaces")], className="home-metric-tile"),
                    html.Div([html.Strong("3", className="home-metric-value"), html.Span("optical system components")], className="home-metric-tile"),
                    html.Div([html.Strong("CSV", className="home-metric-value"), html.Span("structured experiment export")], className="home-metric-tile"),
                ],
            ),
        ],
    )
