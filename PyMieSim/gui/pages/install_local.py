"""Local installation documentation page."""

from dash import html

from PyMieSim.gui.components import Card


RELEASES_URL = "https://github.com/MartinPdeS/PyMieSim/releases"


def _documentation_card(title: str, subtitle: str, body):
    """Build a RosettaX-style documentation card."""
    return html.Section(
        className=Card.classes(color="blue", extra="panel documentation-card"),
        children=[
            html.Div(
                className="card-header panel-header",
                children=[
                    html.Div(title, className="card-title"),
                    html.Div(subtitle, className="card-subtitle"),
                ],
            ),
            html.Div(className="card-body", children=body),
        ],
    )


def _link_chip(label: str, href: str):
    return html.A(label, href=href, className="documentation-link-chip")


def build_install_local_page():
    """Build the local-installation page linked from the home card."""
    return html.Div(
        className="page-content-stack",
        children=[
            html.Section(
                className="page-hero documentation-page-hero",
                children=[
                    html.P("PyMieSim reference", className="eyebrow"),
                    html.H1("Install PyMieSim Locally"),
                    html.P(
                        "This page explains how to download and run PyMieSim from its release distribution without building from source.",
                        className="hero-text",
                    ),
                ],
            ),
            _documentation_card(
                "Related pages",
                "Use the documentation hub for workflow and model details after installation.",
                [
                    html.Div(
                        [_link_chip("Back to documentation hub", "/documentation"), _link_chip("Quick start", "/documentation")],
                        className="documentation-link-chips",
                    )
                ],
            ),
            html.Div(
                className="install-local-grid",
                children=[
                    _documentation_card(
                        "Download latest release bundle",
                        "Choose the archive for your operating system.",
                        [
                            html.Div(
                                [
                                    html.A("Download Windows", href=RELEASES_URL, target="_blank", rel="noopener noreferrer", className="install-download-button"),
                                    html.A("Download macOS", href=RELEASES_URL, target="_blank", rel="noopener noreferrer", className="install-download-button"),
                                    html.A("Download Linux", href=RELEASES_URL, target="_blank", rel="noopener noreferrer", className="install-download-button"),
                                ],
                                className="install-download-buttons",
                            ),
                            html.Div(
                                [
                                    "Need a different version? Use the full ",
                                    html.A("Releases page", href=RELEASES_URL, target="_blank", rel="noopener noreferrer"),
                                    ".",
                                ],
                                className="install-release-note",
                            ),
                        ],
                    ),
                    _documentation_card(
                        "Install and run",
                        "Run PyMieSim after installing the package.",
                        [
                            html.Ol(
                                [
                                    html.Li("Download the archive for your operating system."),
                                    html.Li("Extract it to a folder where you have write permission."),
                                    html.Li("Install PyMieSim in your Python environment."),
                                    html.Li("Launch the PyMieSim GUI and open the local dashboard."),
                                ],
                                className="install-steps",
                            )
                        ],
                    ),
                ],
            ),
            _documentation_card(
                "Validation and provenance",
                "Recommended checks before using a local installation.",
                [
                    html.Div("Confirm you downloaded the expected version from the release page."),
                    html.Div("Keep the installed version together with reports for traceability."),
                    html.Div("For source-based installation, follow the instructions in the main documentation."),
                ],
            ),
        ],
    )
