"""Layout helpers for the experiment dashboard."""

from __future__ import annotations

from dash import dcc, html

from PyMieSim.gui.components import HeaderCard
from PyMieSim.gui.schemas import FieldSpec, SECTION_FIELDS, SINGLE_SCATTERER_FIELDS, SINGLE_SOURCE_FIELDS


def _build_legacy_layout(default_measure_options: list[str]):
    """Build the former all-workspaces layout used to extract page content."""
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
                                value="experiment-tab",
                                className="main-tabs workspace-hidden",
                                children=[
                                    dcc.Tab(
                                        label="Experiment",
                                        value="experiment-tab",
                                        className="main-tab",
                                        selected_className="main-tab--selected",
                                        children=[
                            HeaderCard(
                                "Experiment Lab",
                                "Configure source, scatterer, and detector sets, then run the compiled experiment engine directly from Dash.",
                                [
                                    ("01", "Configure source", "Choose the source family and sweep its optical parameters.", "green"),
                                    ("02", "Configure scatterer", "Set particle geometry, material, and medium values.", "blue"),
                                    ("03", "Configure detector", "Select collection geometry and compute the experiment response.", "yellow"),
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
                                            html.Section(
                                                className="panel plot-controls-panel",
                                                children=[
                                                    html.Div(className="panel-header", children=[html.H2("Plot Controls")]),
                                                    html.Div(
                                                        className="run-controls-grid",
                                                        children=[
                                                            html.Div(
                                                                className="field-block",
                                                                children=[
                                                                    html.Label("Measure", htmlFor="measure-select"),
                                                                    dcc.Dropdown(
                                                                        id="measure-select",
                                                                        className="dashboard-dropdown",
                                                                        options=[{"label": measure, "value": measure} for measure in default_measure_options],
                                                                        value=default_measure_options[0] if default_measure_options else None,
                                                                        clearable=False,
                                                                        optionHeight=34,
                                                                        maxHeight=200,
                                                                    ),
                                                                ],
                                                            ),
                                                            html.Div(
                                                                className="field-block",
                                                                children=[
                                                                    html.Label("X Axis", htmlFor="x-axis-select"),
                                                                    dcc.Dropdown(id="x-axis-select", className="dashboard-dropdown", options=[], placeholder="Detected from fields with multiple values", optionHeight=34, maxHeight=200),
                                                                ],
                                                            ),
                                                        ],
                                                    ),
                                                    html.Div(
                                                        className="run-actions plot-control-actions",
                                                        children=[
                                                            html.Button("Run", id="run-experiment", n_clicks=0, className="run-button run-button-primary"),
                                                            html.Button("Export CSV", id="export-csv", n_clicks=0, className="run-button export-button"),
                                                        ],
                                                    ),
                                                    html.Div(id="status-banner", className="status-banner idle", children="Ready."),
                                                ],
                                            ),
                                        ],
                                    ),
                                    html.Section(
                                        className="result-column",
                                        children=[
                                            html.Section(
                                                className="panel graph-panel",
                                                children=[dcc.Graph(id="result-graph", config={"displaylogo": False})],
                                            ),
                                            html.Section(
                                                id="guide",
                                                className="panel helper-panel",
                                                children=[
                                                    html.Div(className="panel-header", children=[html.H2("Field Syntax Guide")]),
                                                    html.P(
                                                        "Every text field accepts compact batch input, so you can sweep parameters without rewriting the form.",
                                                        className="helper-copy",
                                                    ),
                                                    html.Div(
                                                        className="helper-examples",
                                                        children=[
                                                            html.Div(className="helper-chip", children=[html.Strong("a:b:c"), html.Span("Generate a sweep from start to stop using exactly c points, for example 400:1400:8.")]),
                                                            html.Div(className="helper-chip", children=[html.Strong("1,2,3"), html.Span("Enter a manual list when you want precise sampled values instead of an evenly spaced sweep.")]),
                                                            html.Div(className="helper-chip", children=[html.Strong("LP01,HG11"), html.Span("Mode and label fields also accept comma-separated names for detector families or custom batches.")]),
                                                            html.Div(className="helper-chip helper-chip-note", children=[html.Strong("Tip"), html.Span("Keep the X axis on the parameter you actually swept if you want clean single-panel plots.")]),
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
                                    dcc.Tab(
                                        label="Single / Representations",
                                        value="single-tab",
                                        className="main-tab",
                                        selected_className="main-tab--selected",
                                        children=[_build_single_tab()],
                                    ),
                                ],
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


def build_workspace_layout(default_measure_options: list[str], tab_value: str):
    """Extract only the selected workspace page from the legacy composition."""
    legacy_layout = _build_legacy_layout(default_measure_options)
    tabs = legacy_layout.children[4].children[1].children[2]
    selected_tab = next(tab for tab in tabs.children if tab.value == tab_value)
    return html.Div(className="workspace-page", children=selected_tab.children)


def create_layout(default_measure_options: list[str]):
    """Create a route-driven application shell containing only Home initially."""
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
                        id="page-content",
                        className="dashboard-main",
                        children=[_build_home_page()],
                    ),
                ],
            ),
        ],
    )


def _build_single_tab():
    """Build the single-scatterer representation workspace."""
    from PyMieSim.gui.pages.single import build_single_page

    return build_single_page()

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
                            html.H1("Single Representation Studio"),
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
                                dcc.Dropdown(id="single-source-type", className="dashboard-dropdown", options=[{"label": key, "value": key} for key in SINGLE_SOURCE_FIELDS], value="Gaussian", clearable=False),
                                html.Div(id="single-source-fields"),
                            ),
                            _section_shell(
                                "Scatterer",
                                dcc.Dropdown(id="single-scatterer-type", className="dashboard-dropdown", options=[{"label": key, "value": key} for key in SINGLE_SCATTERER_FIELDS], value="Sphere", clearable=False),
                                html.Div(id="single-scatterer-fields"),
                            ),
                            html.Section(
                                className="panel run-panel",
                                children=[
                                    html.Div(className="panel-header", children=[html.H2("Representation Controls")]),
                                    html.Div(className="field-block", children=[html.Label("Representation", htmlFor="single-representation"), dcc.Dropdown(id="single-representation", className="dashboard-dropdown", options=[{"label": "S1 / S2 amplitudes", "value": "s1s2"}, {"label": "Stokes intensity", "value": "stokes"}, {"label": "Scattering phase function", "value": "spf"}, {"label": "Far-field intensity", "value": "farfields"}], value="s1s2", clearable=False)]),
                                    html.Div(className="field-block", children=[html.Label("Angular sampling", htmlFor="single-sampling"), dcc.Input(id="single-sampling", type="number", value=120, min=24, max=300, step=1, className="field-input")]),
                                    html.Button("Render representation", id="run-single", n_clicks=0, className="run-button run-button-primary"),
                                    html.Div(id="single-status", className="status-banner idle", children="Ready."),
                                ],
                            ),
                        ],
                    ),
                    html.Section(
                        className="result-column",
                        children=[html.Div(id="single-summary", className="summary-grid"), html.Section(className="panel graph-panel single-graph-panel", children=[dcc.Graph(id="single-graph", config={"displaylogo": False})])],
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
            html.Section(className="page-hero", children=[html.P("PyMieSim", className="eyebrow"), html.H1("Optical scattering, clearly organized."), html.P("Explore parameter sweeps in the Experiment workspace or inspect a single particle through its physical representations.", className="hero-text")]),
            html.Div(
                className="home-card-grid",
                children=[
                    html.A(href="/experiment", className="home-action-card", children=[html.Span("01", className="home-card-number"), html.H2("Experiment"), html.P("Run source, scatterer, and detector sets across parameter sweeps. Export structured results for analysis.")]),
                    html.A(href="/single", className="home-action-card", children=[html.Span("02", className="home-card-number"), html.H2("Single particle"), html.P("Render S1/S2, Stokes, SPF, and far-field representations for one optical setup.")]),
                    html.A(href="/documentation", className="home-action-card", children=[html.Span("03", className="home-card-number"), html.H2("Documentation"), html.P("Learn the model vocabulary, field syntax, supported objects, and recommended workflows.")]),
                ],
            ),
            html.Section(className="panel home-overview-card", children=[html.Div(className="panel-header", children=[html.H2("A compact front door to PyMieSim")]), html.P("The dashboard follows the same calm, sectioned layout language as Rosettax: navigation stays in the sidebar, while each workspace gets room for its own controls and visual output."), html.Div(className="meta-strip", children=[html.Div(className="meta-chip", children=[html.Span("Compiled engine"), html.Small("C++ single and experiment backends")]), html.Div(className="meta-chip", children=[html.Span("Plot-ready"), html.Small("Plotly figures and CSV export")])])]),
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
                    _documentation_card("Experiment workflow", "Use source sets, scatterer sets, and detector sets to build sweeps. The X-axis selector is inferred from fields containing multiple values."),
                    _documentation_card("Single workflow", "Use one source and one scatterer to inspect angular amplitudes, polarization, phase functions, or far-field intensity."),
                    _documentation_card("Field syntax", "Quantity inputs accept scalar values, comma-separated lists, and start:end:count expressions such as 400:1400:8."),
                    _documentation_card("Supported models", "Gaussian and plane-wave sources are available in both workflows, alongside spheres, infinite cylinders, and core-shell scatterers."),
                ],
            ),
            html.Section(className="panel documentation-note", children=[html.Div(className="panel-header", children=[html.H2("Where to start")]), html.P("Choose Experiment for parameter studies and detector coupling. Choose Single for representation plots and physical intuition. Both pages preserve the same source/scatterer vocabulary so configurations transfer naturally between them."), html.A("Open the Experiment workspace →", href="/experiment", className="inline-action")]),
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
                        html.Img(src="/assets/pymiesim-logo.svg", alt="PyMieSim home", className="sidebar-logo"),
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
                    _sidebar_link("Experiment", "/experiment"),
                    _sidebar_link("Single", "/single"),
                    _sidebar_link("Documentation", "/documentation"),
                ],
            ),
        ],
    )


def _sidebar_link(label: str, href: str):
    """Build a single sidebar anchor."""
    link_id = f"sidebar-link-{label.lower().replace(' ', '-')}"
    return dcc.Link(label, id=link_id, href=href, refresh=False, className="sidebar-link active" if label == "Home" else "sidebar-link")


def render_fields(section: str, section_type: str):
    """Render all fields for one selected section type."""
    schemas = {"single-source": SINGLE_SOURCE_FIELDS, "single-scatterer": SINGLE_SCATTERER_FIELDS}
    field_specs = schemas[section][section_type] if section in schemas else SECTION_FIELDS[section][section_type]

    if not field_specs:
        return html.Div("This experiment runs without a detector set.", className="empty-section")

    return html.Div(className="field-grid", children=[render_field(section, field_spec) for field_spec in field_specs])


def render_field(section: str, field_spec: FieldSpec):
    """Render one schema field as a labeled text input."""
    return html.Div(
        className="field-block",
        children=[
            html.Label(_format_field_label(field_spec)),
            dcc.Input(
                id={"kind": "field", "section": section, "name": field_spec.name},
                type="text",
                value=field_spec.default,
                debounce=True,
                placeholder=field_spec.placeholder,
                className="field-input",
            ),
        ],
    )


def render_summary_cards(summary: list[dict[str, str]]):
    """Render compact summary cards for one experiment run."""
    if not summary:
        return [html.Div(className="summary-card", children=[html.Span("No result yet"), html.Strong("Run an experiment")])]

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
        return "Single values, comma-separated values, or start:end:count are supported when meaningful."

    return f"Values are interpreted in {field_spec.unit}."


def _format_field_label(field_spec: FieldSpec) -> str:
    """Put the unit beside the label so cards stay compact."""
    if field_spec.unit is None:
        return field_spec.label
    return f"{field_spec.label} [{field_spec.unit}]"
