"""Layout helpers for the experiment dashboard."""

from __future__ import annotations

from dash import dcc, html

from PyMieSim.gui.components import HeaderCard
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


def _build_legacy_layout(default_measure_options: list[str], plot_settings: dict | None = None):
    """Build the former all-workspaces layout used to extract page content."""
    stored = plot_settings or {}
    particle_settings = stored.get("particle_explorer", {})
    sweep_settings = stored.get("parameter_sweep", {})
    return html.Div(
        className="app-shell",
        children=[
            dcc.Location(id="url", refresh=False),
            dcc.Store(id="experiment-result"),
            dcc.Store(id="single-result"),
            dcc.Download(id="csv-download"),
            html.Div(
                className="dashboard-frame",
                children=[
                    _build_sidebar(),
                    html.Main(
                        className="dashboard-main",
                        children=[
                            html.Div(id="home-page", className="page-view", style={"display": "block"}, children=[_build_home_page()]),
                            html.Div(id="documentation-page", className="page-view", style={"display": "none"}, children=[_build_documentation_page()]),
                            dcc.Tabs(
                                id="main-tabs",
                                value=DEFAULT_APPLICATION_SETTINGS["default_tab"],
                                className="main-tabs workspace-hidden",
                                children=[
                                    dcc.Tab(
                                        label="Particle Explorer / Representations",
                                        value="single-tab",
                                        className="main-tab",
                                        selected_className="main-tab--selected",
                                        children=[_build_single_tab(particle_settings)],
                                    ),
                                    dcc.Tab(
                                        label="Parameter Sweep",
                                        value="experiment-tab",
                                        className="main-tab",
                                        selected_className="main-tab--selected",
                                        children=[
                            HeaderCard(
                                "Parameter Sweep Lab",
                                "Configure source, scatterer, and detector sets, then run parameter sweeps through the compiled engine directly from Dash.",
                                [
                                    ("01", "Configure source", "Choose the source family and sweep its optical parameters.", "yellow"),
                                    ("02", "Configure scatterer", "Set particle geometry, material, and medium values.", "blue"),
                                    ("03", "Configure detector", "Select collection geometry and compute the sweep response.", "orange"),
                                ],
                                color="green",
                            ).render(),
                            html.Section(
                                id="configure",
                                className="workspace",
                                children=[
                                    html.Section(
                                        className="control-column",
                                        children=[
                                            html.Div(
                                                className="set-panel-grid",
                                                children=[
                                                    *_build_experiment_configuration_sections(),
                                                ],
                                            ),
                                        ],
                                    ),
                                    html.Section(
                                        className="result-column",
                                        children=[
                                            html.Section(
                                                className="panel graph-panel",
                                                children=[dcc.Graph(id="result-graph", config=PLOT_CONFIG)],
                                            ),
                                            _x_axis_card(default_measure_options),
                                            html.Div(id="experiment-plot-options-container", children=[_plot_options_card("experiment", sweep_settings)]),
                                            html.Div(
                                                className="export-actions",
                                                children=[html.Button("Export CSV", id="export-csv", n_clicks=0, className="run-button export-button")],
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


def build_workspace_layout(default_measure_options: list[str], tab_value: str, plot_settings: dict | None = None):
    """Extract only the selected workspace page from the legacy composition."""
    legacy_layout = _build_legacy_layout(default_measure_options, plot_settings)
    dashboard_frame = next(child for child in legacy_layout.children if getattr(child, "className", None) == "dashboard-frame")
    tabs = dashboard_frame.children[1].children[2]
    selected_tab = next(tab for tab in tabs.children if tab.value == tab_value)
    return html.Div(className="workspace-page", children=selected_tab.children)


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


def _build_single_tab(plot_settings: dict | None = None):
    """Build the single-scatterer representation workspace."""
    from PyMieSim.gui.pages.single import build_single_page

    return build_single_page(plot_settings=plot_settings)

    # Kept below as a reference while the page migration is completed.
    return html.Div(
        className="tab-content-stack single-tab-content",
        children=[
            html.Section(
                className="workflow-page-header hero single-hero",
                children=[
                    html.Div(
                        className="hero-copy",
                        children=[
                            html.P("PyMieSim", className="eyebrow"),
                            html.H1("Particle Explorer"),
                            html.P("Inspect angular scattering, polarization, and field patterns for one source–scatterer setup.", className="hero-text"),
                        ],
                    ),
                ],
            ),
            html.Section(
                className="single-workspace",
                children=[
                    html.Section(
                        className="control-column",
                        children=[
                            _section_shell(
                                "Source",
                                dcc.Dropdown(id="single-source-type", className="dashboard-dropdown", options=[{"label": key, "value": key} for key in SINGLE_SOURCE_FIELDS], value="Gaussian", clearable=False, searchable=False, optionHeight=38, maxHeight=200, persistence=True, persistence_type="session"),
                                html.Div(id="single-source-fields"),
                            ),
                            _section_shell(
                                "Scatterer",
                                dcc.Dropdown(id="single-scatterer-type", className="dashboard-dropdown", options=[{"label": key, "value": key} for key in SINGLE_SCATTERER_FIELDS], value="Sphere", clearable=False, searchable=False, optionHeight=38, maxHeight=200, persistence=True, persistence_type="session"),
                                html.Div(id="single-scatterer-fields"),
                            ),
                            html.Section(
                                className="panel run-panel",
                                children=[
                                    html.Div(className="panel-header", children=[html.H2("Representation Controls")]),
                                    html.Div(className="field-block", children=[html.Label("Representation", htmlFor="single-representation"), dcc.Dropdown(id="single-representation", className="dashboard-dropdown", options=[{"label": "S1 / S2 amplitudes", "value": "s1s2"}, {"label": "Stokes I intensity", "value": "stokes"}, {"label": "Stokes Q", "value": "stokes_q"}, {"label": "Stokes U", "value": "stokes_u"}, {"label": "Stokes V", "value": "stokes_v"}, {"label": "Scattering phase function", "value": "spf"}, {"label": "Far-field intensity", "value": "farfields"}], value="s1s2", clearable=False, searchable=False, optionHeight=38, maxHeight=200, persistence=True, persistence_type="session")]),
                                    html.Div(className="field-block", children=[html.Label("Angular sampling", htmlFor="single-sampling"), dcc.Input(id="single-sampling", type="number", value=120, min=24, max=300, step=1, placeholder="120", className="field-input", persistence=True, persistence_type="session")]),
                                ],
                            ),
                        ],
                    ),
                    html.Section(
                        className="result-column",
                        children=[html.Div(id="single-summary", className="summary-grid"), html.Section(className="panel graph-panel single-graph-panel", children=[dcc.Graph(id="single-graph", config=PLOT_CONFIG)]), _plot_options_card("single", particle_settings)],
                    ),
                ],
            ),
        ],
    )


def _build_home_page():
    """Build the Rosettax-style landing page."""
    from PyMieSim.gui.pages.home import build_home_page

    return build_home_page()

    # Legacy inline composition retained temporarily for compatibility.
    return html.Div(
        className="page-content-stack",
        children=[
            html.Section(className="page-hero", children=[html.P("PyMieSim", className="eyebrow"), html.H1("Optical scattering, clearly organized."), html.P("Run parameter sweeps in the Parameter Sweep workspace or inspect one particle through its physical representations in Particle Explorer.", className="hero-text")]),
            html.Div(
                className="home-card-grid",
                children=[
                    html.A(href="/experiment", className="home-action-card", children=[html.Span("01", className="home-card-number"), html.H2("Parameter Sweep"), html.P("Run source, scatterer, and detector sets across parameter sweeps. Export structured results for analysis.")]),
                    html.A(href="/single", className="home-action-card", children=[html.Span("02", className="home-card-number"), html.H2("Particle Explorer"), html.P("Render S1/S2, Stokes, SPF, and far-field representations for one optical setup.")]),
                    html.A(href="/documentation", className="home-action-card", children=[html.Span("03", className="home-card-number"), html.H2("Documentation"), html.P("Learn the model vocabulary, field syntax, supported objects, and recommended workflows.")]),
                ],
            ),
            html.Section(className="panel home-overview-card", children=[html.Div(className="panel-header", children=[html.H2("A compact front door to PyMieSim")]), html.P("The dashboard follows the same calm, sectioned layout language as Rosettax: navigation stays in the sidebar, while Parameter Sweep and Particle Explorer each get room for their own controls and visual output."), html.Div(className="meta-strip", children=[html.Div(className="meta-chip", children=[html.Span("Compiled engine"), html.Small("C++ particle and sweep backends")]), html.Div(className="meta-chip", children=[html.Span("Plot-ready"), html.Small("Plotly figures and CSV export")])])]),
        ],
    )


def _build_documentation_page():
    """Build a lightweight documentation hub matching Rosettax's page style."""
    from PyMieSim.gui.pages.documentation import build_documentation_page

    return build_documentation_page()

    # Legacy inline composition retained temporarily for compatibility.
    return html.Div(
        className="page-content-stack",
        children=[
            html.Section(className="page-hero", children=[html.P("PyMieSim reference", className="eyebrow"), html.H1("Documentation"), html.P("A practical map of the objects, workflows, and input conventions used throughout the dashboard.", className="hero-text")]),
            html.Div(
                className="documentation-grid",
                children=[
                    _documentation_card("Parameter Sweep workflow", "Use source sets, scatterer sets, and detector sets to build sweeps. The X-axis selector is inferred from fields containing multiple values."),
                    _documentation_card("Particle Explorer workflow", "Use one source and one scatterer to inspect angular amplitudes, polarization, phase functions, or far-field intensity."),
                    _documentation_card("Field syntax", "Quantity inputs accept scalar values, comma-separated lists, and start:end:count expressions such as 400:1400:8."),
                    _documentation_card("Supported models", "Gaussian and plane-wave sources are available in both workflows, alongside spheres, infinite cylinders, and core-shell scatterers."),
                ],
            ),
            html.Section(className="panel documentation-note", children=[html.Div(className="panel-header", children=[html.H2("Where to start")]), html.P("Choose Parameter Sweep for parameter studies and detector coupling. Choose Particle Explorer for representation plots and physical intuition. Both pages preserve the same source/scatterer vocabulary so configurations transfer naturally between them."), html.A("Open the Parameter Sweep workspace →", href="/experiment", className="inline-action")]),
        ],
    )


def _documentation_card(title: str, description: str):
    return html.Section(className="panel documentation-card", children=[html.Div(className="panel-header", children=[html.H2(title)]), html.P(description)])


def _build_experiment_configuration_sections():
    """Compose experiment configuration cards from page section modules."""
    from PyMieSim.gui.pages.experiment.sections import (
        build_detector_section,
        build_scatterer_section,
        build_source_section,
    )

    return [build_source_section(), build_scatterer_section(), build_detector_section()]


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


def _section_shell(title: str, selector, content):
    """Wrap one configuration section in consistent panel markup."""
    return html.Section(
        className="panel",
        children=[
            html.Div(className="panel-header", children=[html.H2(title)]),
            selector,
            html.Div(className="panel-body", children=[content]),
        ],
    )


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
