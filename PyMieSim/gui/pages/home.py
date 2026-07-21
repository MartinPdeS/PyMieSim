"""Home page composition."""

from importlib.metadata import PackageNotFoundError, version

from dash import html

from PyMieSim.gui.components import Card

def build_home_page(metrics: dict[str, int] | None = None):
    """Build the landing page and its workflow cards."""
    metrics = metrics or {}
    return html.Div(
        className="page-content-stack",
        children=[
            html.Div(
                [
                    html.Span("Version:", className="version-widget-label"),
                    html.Span(_package_version(), className="version-widget-value"),
                ],
                className="version-widget",
            ),
            html.Section(
                className="page-hero home-page-hero",
                children=[
                    html.H1("PyMieSim"),
                    html.P(
                        "An open-source library for fast and flexible far-field Mie scattering simulations.",
                        className="hero-text",
                    ),
                ],
            ),
            _citation_card(),
            html.Div(
                className="home-capability-grid",
                children=[
                    _capability_card(
                        "Particle Explorer",
                        "Inspect one optical setup through angular scattering, polarization, phase functions, and field representations.",
                        ["Configure source", "Configure scatterer", "Render representations"],
                        "Open Particle Explorer",
                        "/single",
                    ),
                    _capability_card(
                        "Parameter Sweep",
                        "Run source, scatterer, and detector sets across parameter sweeps and export structured results for analysis.",
                        ["Configure source", "Configure scatterer and detector", "Run and export results"],
                        "Open Parameter Sweep",
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
            _metrics_card(metrics),
        ],
    )


def _package_version() -> str:
    try:
        return version("PyMieSim")
    except PackageNotFoundError:
        return "development"


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
                            html.A("Support Developer", href="https://github.com/sponsors/MartinPdeS", target="_blank", rel="noopener noreferrer", className="home-button home-button-primary"),
                            html.A("Citing this work", href="/citation", className="home-button home-button-outline"),
                            html.A("PyMieSim documentation", href="https://martinpdes.github.io/PyMieSim/", target="_blank", rel="noopener noreferrer", className="home-button home-button-info"),
                            html.A("GitHub repository", href="https://github.com/MartinPdeS/PyMieSim", target="_blank", rel="noopener noreferrer", className="home-button home-button-muted"),
                            html.A("Install locally (Releases)", href="/documentation/install-local", className="home-button home-button-muted"),
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


def _metrics_card(metrics: dict[str, int]):
    return html.Section(
        className=Card.classes(color="blue", extra="home-info-card home-metrics-card"),
        children=[
            html.Div(className="card-header panel-header", children=[html.Div("PyMieSim usage metrics.", className="card-title")]),
            html.Div(
                className="card-body",
                children=[
                    html.Div(
                        [
                            html.Div(str(metrics.get("home_page_visits", 0)), className="home-metric-value"),
                            html.Div("Home page visits", className="home-metric-label"),
                        ],
                        className="home-metric-tile",
                    ),
                    html.Div(
                        [
                            html.Div(str(metrics.get("experiment_runs", 0)), className="home-metric-value"),
                            html.Div("Parameter sweeps", className="home-metric-label"),
                        ],
                        className="home-metric-tile",
                    ),
                    html.Div(
                        [
                            html.Div(str(metrics.get("single_runs", 0)), className="home-metric-value"),
                            html.Div("Particle explorations", className="home-metric-label"),
                        ],
                        className="home-metric-tile",
                    ),
                ],
            ),
        ],
    )
