"""Layout helpers for the experiment dashboard."""

from __future__ import annotations

from dash import dcc, html

from PyMieSim.gui.defaults import DEFAULT_APPLICATION_SETTINGS, DEFAULT_PLOT_SETTINGS
from PyMieSim.gui.schemas import FieldSpec, SECTION_FIELDS, SINGLE_SCATTERER_FIELDS, SINGLE_SOURCE_FIELDS


THEME_LIGHT = "https://cdn.jsdelivr.net/npm/bootswatch@5.3.6/dist/flatly/bootstrap.min.css"
THEME_DARK = "https://cdn.jsdelivr.net/npm/bootswatch@5.3.6/dist/slate/bootstrap.min.css"

PLOT_CONFIG = {
    "displayModeBar": "always",
    "displaylogo": False,
    "scrollZoom": True,
    "doubleClick": "reset+autosize",
    "modeBarButtonsToAdd": ["zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d"],
    "toImageButtonOptions": {"format": "png", "filename": "pymiesim_plot", "scale": 2},
}


def _plot_options_card(prefix: str, settings: dict, representation: str | None = None, projection: str | None = None, projection_options: list[dict] | None = None):
    """Build the neutral, per-plot controls shown beneath each graph."""
    values = {
        "x_scale": "linear",
        "log_y": False,
        "font_size": 14,
        "line_width": 2,
        "show_legend": True,
        "show_grid": True,
        "show_title": True,
        **(settings or {}),
    }

    def dropdown(name: str, options: list[dict], value):
        return dcc.Dropdown(
            id=f"plot-{prefix}-{name}",
            options=options,
            value=value,
            clearable=False,
            searchable=False,
            optionHeight=34,
            maxHeight=170,
            persistence=True,
            persistence_type="session",
            className="dashboard-dropdown plot-option-control",
        )

    def number(name: str, value, minimum, maximum, step):
        return dcc.Input(
            id=f"plot-{prefix}-{name}",
            type="number",
            value=value,
            min=minimum,
            max=maximum,
            step=step,
            placeholder=str(value),
            persistence=True,
            persistence_type="session",
            className="field-input plot-option-control plot-number-control",
        )

    fields = []
    is_structured_map = representation in {"stokes", "stokes_q", "stokes_u", "stokes_v", "spf", "farfields"}
    is_polar = prefix == "single" and projection == "polar_1d"
    is_3d = prefix == "single" and projection in {"3d", "3d_radial"}

    if is_polar:
        fields.append(("Radial scale", dropdown("y-scale", [{"label": "Linear", "value": False}, {"label": "Logarithmic", "value": True}], bool(values["log_y"]))))
    elif not is_structured_map and not is_3d:
        fields.extend([
            ("X scale", dropdown("x-scale", [{"label": "Linear", "value": "linear"}, {"label": "Logarithmic", "value": "log"}], values["x_scale"])),
            ("Y scale", dropdown("y-scale", [{"label": "Linear", "value": False}, {"label": "Logarithmic", "value": True}], bool(values["log_y"]))),
        ])
    fields.append(("Font size", number("font-size", values["font_size"], 8, 32, 1)))
    if not is_structured_map and not is_3d:
        fields.append(("Line width", number("line-width", values["line_width"], 0.5, 8, 0.5)))
    fields.extend([
        ("Legend", dropdown("legend", [{"label": "Show", "value": True}, {"label": "Hide", "value": False}], bool(values["show_legend"]))),
        ("Grid", dropdown("grid", [{"label": "Show", "value": True}, {"label": "Hide", "value": False}], bool(values["show_grid"]))),
    ])
    if prefix == "experiment":
        fields.insert(2, ("Projection", dropdown("projection", projection_options or [{"label": "Cartesian", "value": "cartesian"}], projection or "cartesian")))
    return html.Section(
        className="plot-options-card",
        children=[
            html.Div(className="plot-options-title", children="Plot options"),
            html.Div(className="plot-options-grid", children=[html.Div(className="plot-option-field", children=[html.Label(label), control]) for label, control in fields]),
        ],
    )


def _x_axis_card(default_measure_options: list[str]):
    """Build the post-plot X/Y axis selectors, separate from run controls."""
    return html.Section(
        className="graph-axis-card",
        children=[
            html.Div("Graph axes", className="graph-axis-title"),
            html.Div(
                className="graph-axis-fields",
                children=[
                    html.Div(
                        className="graph-axis-field",
                        children=[
                            html.Label("X axis", htmlFor="x-axis-select"),
                            dcc.Dropdown(
                                id="x-axis-select",
                                className="dashboard-dropdown graph-axis-control",
                                options=[],
                                placeholder="Detected from fields with multiple values",
                                searchable=False,
                                optionHeight=38,
                                maxHeight=200,
                                persistence=True,
                                persistence_type="session",
                            ),
                        ],
                    ),
                    html.Div(
                        className="graph-axis-field",
                        children=[
                            html.Label("Y axis", htmlFor="measure-select"),
                            dcc.Dropdown(
                                id="measure-select",
                                className="dashboard-dropdown graph-axis-control",
                                options=[{"label": measure, "value": measure} for measure in default_measure_options],
                                value=default_measure_options[0] if default_measure_options else None,
                                clearable=False,
                                searchable=False,
                                optionHeight=38,
                                maxHeight=200,
                                persistence=True,
                                persistence_type="session",
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


def create_layout(default_measure_options: list[str]):
    """Create a route-driven application shell containing only Home initially."""
    return html.Div(
        className="app-shell",
        children=[
            dcc.Location(id="url", refresh=False),
            dcc.Store(id="experiment-result"),
            dcc.Store(id="single-result"),
            dcc.Download(id="csv-download"),
            dcc.Download(id="single-csv-download"),
            dcc.Store(id="home-visit-count", data=0, storage_type="local"),
            dcc.Store(id="experiment-run-count", data=0, storage_type="local"),
            dcc.Store(id="single-run-count", data=0, storage_type="local"),
            dcc.Store(id="theme-store", data={"theme": DEFAULT_APPLICATION_SETTINGS["theme"]}, storage_type="session"),
            dcc.Store(id="plot-settings-store", data={
                "particle_explorer": dict(DEFAULT_PLOT_SETTINGS),
                "parameter_sweep": dict(DEFAULT_PLOT_SETTINGS),
            }, storage_type="session"),
            html.Link(id="theme-link", rel="stylesheet", href=THEME_LIGHT),
            html.Div(
                className="dashboard-frame",
                children=[
                    _build_sidebar(),
                    html.Main(
                        id="page-content",
                        className="dashboard-main",
                        children=[],
                    ),
                ],
            ),
        ],
    )


def build_application_footer() -> html.Footer:
    """Build the global attribution footer shown below every page."""
    return html.Footer(
        [
            "RosettaX is developed and maintained by ",
            html.A(
                "Martin Poinsinet de Sivry-Houle",
                href="mailto:martin.poinsinet.de.sivry@gmail.com",
                style={
                    "textDecoration": "none",
                    "fontWeight": "600",
                },
            ),
            ", at the Amsterdam Vesicle Center.",
        ],
        style={
            "fontSize": "0.88rem",
            "opacity": 0.72,
            "textAlign": "center",
            "paddingTop": "12px",
            "paddingBottom": "8px",
        },
    )


def build_page_with_footer(page):
    """Wrap a routed page with the shared RosettaX-style footer."""
    return html.Div(
        children=[page, build_application_footer()],
        style={
            "flex": "1 0 auto",
            "display": "flex",
            "flexDirection": "column",
            "gap": "24px",
        },
    )


def _build_sidebar():
    """Create the dashboard sidebar."""
    return html.Aside(
        className="dashboard-sidebar",
        children=[
            html.Div(
                className="sidebar-brand",
                children=[
                    dcc.Link(
                        html.Img(id="sidebar-logo", src="/assets/pymiesim-logo.svg", alt="PyMieSim home", className="sidebar-logo"),
                        href="/",
                        refresh=False,
                        className="sidebar-logo-link",
                    ),
                ],
            ),
            html.Nav(
                className="sidebar-nav",
                children=[
                    _sidebar_link("Home", "/"),
                    _sidebar_link("Particle Explorer", "/single"),
                    _sidebar_link("Parameter Sweep", "/experiment"),
                    _sidebar_link("Documentation", "/documentation"),
                    _sidebar_link("Settings", "/settings"),
                ],
            ),
        ],
    )


def build_theme_selector():
    """Build the persisted Light/Dark selector used by the application shell."""
    return dcc.Dropdown(
        id="theme-mode",
        options=[
            {"label": "Dark", "value": "dark"},
            {"label": "Light", "value": "light"},
        ],
        value=DEFAULT_APPLICATION_SETTINGS["theme"],
        clearable=False,
        searchable=False,
        optionHeight=38,
        maxHeight=200,
        persistence=True,
        persistence_type="session",
        className="theme-mode-select",
    )


def _sidebar_link(label: str, href: str):
    """Build a single sidebar anchor."""
    # Keep route-based IDs stable while allowing user-facing labels to evolve.
    link_id = f"sidebar-link-{href.strip('/') or 'home'}"
    return dcc.Link(label, id=link_id, href=href, refresh=False, className="sidebar-link active" if label == "Home" else "sidebar-link")


def render_fields(section: str, section_type: str):
    """Render all fields for one selected section type."""
    schemas = {"single-source": SINGLE_SOURCE_FIELDS, "single-scatterer": SINGLE_SCATTERER_FIELDS}
    field_specs = schemas[section][section_type] if section in schemas else SECTION_FIELDS[section][section_type]

    if not field_specs:
        return html.Div("This parameter sweep runs without a detector set.", className="empty-section")

    return html.Div(className="field-grid", children=[render_field(section, field_spec) for field_spec in field_specs])


def render_field(section: str, field_spec: FieldSpec):
    """Render one schema field as a labeled text input."""
    persistence = "parameter-sweep-defaults-v2" if section in {"source", "scatterer", "detector"} else True
    default_class = " field-input-default" if section == "detector" and field_spec.name in {"polarization_filter", "medium"} else ""
    return html.Div(
        className="field-block",
        children=[
            html.Label(_format_field_label(field_spec)),
            dcc.Input(
                id={"kind": "field", "section": section, "name": field_spec.name},
                type="text",
                value=field_spec.default,
                debounce=False,
                placeholder=field_spec.placeholder or field_spec.default,
                className=f"field-input{default_class}",
                persistence=persistence,
                persistence_type="session",
            ),
        ],
    )


def render_summary_cards(summary: list[dict[str, str]]):
    """Render compact summary cards for one experiment run."""
    if not summary:
        return [html.Div(className="summary-card", children=[html.Span("No result yet"), html.Strong("Run a parameter sweep")])]

    return [
        html.Div(className="summary-card", children=[html.Span(item["label"]), html.Strong(item["value"])])
        for item in summary
    ]


def _build_default_help_text(field_spec: FieldSpec) -> str:
    """Return a fallback field description."""
    if field_spec.unit is None:
        return "Scalar values, comma-separated values, or start:end:count are supported when meaningful."

    return f"Values are interpreted in {field_spec.unit}."


def _format_field_label(field_spec: FieldSpec) -> str:
    """Put the unit beside the label so cards stay compact."""
    if field_spec.unit is None:
        return field_spec.label
    return f"{field_spec.label} [{field_spec.unit}]"
