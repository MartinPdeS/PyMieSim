"""Settings page composition and independent plotting preferences."""

from __future__ import annotations

from dash import dcc, html

from PyMieSim.gui.components import Card


_COMMON_PLOT_SETTINGS = {
    "font_size": 14,
    "line_width": 2,
    "marker_size": 6,
    "template": "match-theme",
    "show_legend": True,
    "show_grid": True,
    "coordinate_system": "cartesian",
    "show_title": True,
    "log_y": False,
}
DEFAULT_PARTICLE_PLOT_SETTINGS = {**_COMMON_PLOT_SETTINGS}
DEFAULT_SWEEP_PLOT_SETTINGS = {**_COMMON_PLOT_SETTINGS}
# Kept as a compatibility alias for callers that used the original global defaults.
DEFAULT_PLOT_SETTINGS = DEFAULT_PARTICLE_PLOT_SETTINGS


def build_settings_page(theme: str = "light", plot_settings: dict | None = None):
    """Build the settings page with separate controls for each workspace."""
    stored = plot_settings or {}
    particle = _workspace_settings(stored, "particle_explorer", DEFAULT_PARTICLE_PLOT_SETTINGS)
    sweep = _workspace_settings(stored, "parameter_sweep", DEFAULT_SWEEP_PLOT_SETTINGS)
    return html.Div(
        className="page-content-stack settings-page",
        children=[
            html.Section(
                className="page-hero documentation-page-hero",
                children=[
                    html.P("Application preferences", className="eyebrow"),
                    html.H1("Settings"),
                    html.P("Tune the dashboard appearance and configure each plotting workspace independently. Changes are saved automatically in this browser.", className="hero-text"),
                ],
            ),
            html.Div(className="settings-grid", children=[_appearance_card(theme), _plot_card("particle", "Particle Explorer", "blue", particle), _plot_card("sweep", "Parameter Sweep", "yellow", sweep)]),
            html.Section(
                className=Card.classes(color="blue", extra="panel settings-note"),
                children=[
                    html.Div(className="card-header panel-header", children=[html.H2("How settings are applied")]),
                    html.Div(className="card-body", children=[html.P("Particle Explorer and Parameter Sweep remember separate plotting preferences. Polar coordinates are applied to trace-based plots when the representation supports them; heatmaps remain Cartesian."), html.Div(id="settings-save-status", className="settings-save-status", children="All settings are saved automatically.")]),
                ],
            ),
        ],
    )


def _workspace_settings(stored: dict, key: str, defaults: dict) -> dict:
    """Read nested preferences while tolerating the previous flat format."""
    values = stored.get(key, stored if key == "particle_explorer" else {})
    return {**defaults, **(values or {})}


def _appearance_card(theme: str):
    return html.Section(
        className=Card.classes(color="blue", extra="panel settings-card"),
        children=[
            html.Div(className="card-header panel-header", children=[html.H2("Appearance")]),
            html.Div(className="card-body settings-card-body", children=[html.Label("Theme", htmlFor="settings-theme-mode", className="settings-label"), _dropdown("settings-theme-mode", [{"label": "Light", "value": "light"}, {"label": "Dark", "value": "dark"}], theme if theme in {"light", "dark"} else "light"), html.P("The theme changes the dashboard shell and its controls.", className="settings-help")]),
        ],
    )


def _plot_card(prefix: str, title: str, color: str, settings: dict):
    return html.Section(
        className=Card.classes(color=color, extra=f"panel settings-card settings-{prefix}-card"),
        children=[
            html.Div(className="card-header panel-header", children=[html.H2(title), html.Span("Independent plotting preferences", className="settings-card-subtitle")]),
            html.Div(
                className="card-body settings-card-body settings-form-grid",
                children=[
                    _number_setting("Font size", f"settings-{prefix}-font-size", settings["font_size"], 8, 32, 1),
                    _number_setting("Line width", f"settings-{prefix}-line-width", settings["line_width"], 0.5, 8, 0.5),
                    _number_setting("Marker size", f"settings-{prefix}-marker-size", settings["marker_size"], 0, 24, 1),
                    _dropdown_field("Plot theme", f"settings-{prefix}-template", [{"label": "Match application theme", "value": "match-theme"}, {"label": "Light plot", "value": "plotly_white"}, {"label": "Dark plot", "value": "plotly_dark"}], settings["template"]),
                    _dropdown_field("Coordinates", f"settings-{prefix}-coordinates", [{"label": "Cartesian", "value": "cartesian"}, {"label": "Polar (when supported)", "value": "polar"}], settings["coordinate_system"]),
                    _dropdown_field("Legend", f"settings-{prefix}-legend", [{"label": "Show", "value": True}, {"label": "Hide", "value": False}], bool(settings["show_legend"])),
                    _dropdown_field("Grid", f"settings-{prefix}-grid", [{"label": "Show", "value": True}, {"label": "Hide", "value": False}], bool(settings["show_grid"])),
                    _dropdown_field("Title", f"settings-{prefix}-title", [{"label": "Show", "value": True}, {"label": "Hide", "value": False}], bool(settings["show_title"])),
                    _dropdown_field("Y-axis scale", f"settings-{prefix}-log-y", [{"label": "Linear", "value": False}, {"label": "Logarithmic", "value": True}], bool(settings["log_y"])),
                ],
            ),
        ],
    )


def _dropdown_field(label: str, component_id: str, options: list[dict], value):
    return html.Div(className="settings-field", children=[html.Label(label, htmlFor=component_id, className="settings-label"), _dropdown(component_id, options, value)])


def _dropdown(component_id: str, options: list[dict], value):
    return dcc.Dropdown(id=component_id, options=options, value=value, clearable=False, optionHeight=38, maxHeight=200, persistence=True, persistence_type="session", className="dashboard-dropdown")


def _number_setting(label: str, component_id: str, value: int | float, minimum: int | float, maximum: int | float, step: int | float):
    return html.Div(className="settings-field", children=[html.Label(label, htmlFor=component_id, className="settings-label"), dcc.Input(id=component_id, type="number", value=value, min=minimum, max=maximum, step=step, persistence=True, persistence_type="session", className="settings-number-input")])
