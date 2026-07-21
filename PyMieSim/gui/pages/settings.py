"""Settings page composition and plotting preference defaults."""

from __future__ import annotations

from dash import dcc, html

from PyMieSim.gui.components import Card


DEFAULT_PLOT_SETTINGS = {
    "font_size": 14,
    "line_width": 2,
    "marker_size": 6,
    "template": "match-theme",
    "show_legend": True,
    "show_grid": True,
}


def build_settings_page(theme: str = "light", plot_settings: dict | None = None):
    """Build the card-based application settings page."""
    settings = {**DEFAULT_PLOT_SETTINGS, **(plot_settings or {})}
    return html.Div(
        className="page-content-stack settings-page",
        children=[
            html.Section(
                className="page-hero documentation-page-hero",
                children=[
                    html.P("Application preferences", className="eyebrow"),
                    html.H1("Settings"),
                    html.P("Tune the dashboard appearance and the way plots are rendered. Changes are saved automatically in this browser.", className="hero-text"),
                ],
            ),
            html.Div(
                className="settings-grid",
                children=[
                    _appearance_card(theme),
                    _plot_style_card(settings),
                    _plot_layout_card(settings),
                ],
            ),
            html.Section(
                className=Card.classes(color="blue", extra="panel settings-note"),
                children=[
                    html.Div(className="card-header panel-header", children=[html.H2("How settings are applied")]),
                    html.Div(className="card-body", children=[html.P("Plot preferences apply the next time an experiment or single representation is rendered. Existing figures update immediately when a preference changes."), html.Div(id="settings-save-status", className="settings-save-status", children="All settings are saved automatically.")]),
                ],
            ),
        ],
    )


def _appearance_card(theme: str):
    return html.Section(
        className=Card.classes(color="blue", extra="panel settings-card"),
        children=[
            html.Div(className="card-header panel-header", children=[html.H2("Appearance")]),
            html.Div(
                className="card-body settings-card-body",
                children=[
                    html.Label("Theme", htmlFor="settings-theme-mode", className="settings-label"),
                    dcc.Dropdown(
                        id="settings-theme-mode",
                        options=[{"label": "Light", "value": "light"}, {"label": "Dark", "value": "dark"}],
                        value=theme if theme in {"light", "dark"} else "light",
                        clearable=False,
                        optionHeight=38,
                        maxHeight=200,
                        className="dashboard-dropdown",
                    ),
                    html.P("The theme changes the dashboard shell and its controls.", className="settings-help"),
                ],
            ),
        ],
    )


def _plot_style_card(settings: dict):
    return html.Section(
        className=Card.classes(color="green", extra="panel settings-card"),
        children=[
            html.Div(className="card-header panel-header", children=[html.H2("Plot style")]),
            html.Div(
                className="card-body settings-card-body settings-form-grid",
                children=[
                    _number_setting("Font size", "settings-font-size", settings["font_size"], 8, 32, 1),
                    _number_setting("Line width", "settings-line-width", settings["line_width"], 0.5, 8, 0.5),
                    _number_setting("Marker size", "settings-marker-size", settings["marker_size"], 0, 24, 1),
                    html.Div(
                        className="settings-field settings-field-wide",
                        children=[
                            html.Label("Plot theme", htmlFor="settings-plot-template", className="settings-label"),
                            dcc.Dropdown(
                                id="settings-plot-template",
                                options=[
                                    {"label": "Match application theme", "value": "match-theme"},
                                    {"label": "Light plot", "value": "plotly_white"},
                                    {"label": "Dark plot", "value": "plotly_dark"},
                                ],
                                value=settings["template"],
                                clearable=False,
                                optionHeight=38,
                                maxHeight=200,
                                className="dashboard-dropdown",
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


def _plot_layout_card(settings: dict):
    return html.Section(
        className=Card.classes(color="yellow", extra="panel settings-card"),
        children=[
            html.Div(className="card-header panel-header", children=[html.H2("Plot layout")]),
            html.Div(
                className="card-body settings-card-body",
                children=[
                    html.Div(
                        className="settings-form-grid settings-layout-fields",
                        children=[
                            html.Div(
                                className="settings-field",
                                children=[
                                    html.Label("Legend", htmlFor="settings-show-legend", className="settings-label"),
                                    dcc.Dropdown(id="settings-show-legend", options=[{"label": "Show", "value": True}, {"label": "Hide", "value": False}], value=bool(settings["show_legend"]), clearable=False, optionHeight=38, maxHeight=200, className="dashboard-dropdown"),
                                ],
                            ),
                            html.Div(
                                className="settings-field",
                                children=[
                                    html.Label("Grid", htmlFor="settings-show-grid", className="settings-label"),
                                    dcc.Dropdown(id="settings-show-grid", options=[{"label": "Show", "value": True}, {"label": "Hide", "value": False}], value=bool(settings["show_grid"]), clearable=False, optionHeight=38, maxHeight=200, className="dashboard-dropdown"),
                                ],
                            ),
                        ],
                    ),
                    html.P("These options affect both experiment plots and single-particle representations.", className="settings-help"),
                ],
            ),
        ],
    )


def _number_setting(label: str, component_id: str, value: int | float, minimum: int | float, maximum: int | float, step: int | float):
    return html.Div(
        className="settings-field",
        children=[
            html.Label(label, htmlFor=component_id, className="settings-label"),
            dcc.Input(id=component_id, type="number", value=value, min=minimum, max=maximum, step=step, className="settings-number-input"),
        ],
    )
