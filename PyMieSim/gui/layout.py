"""Layout helpers for the experiment dashboard."""

from __future__ import annotations

from dash import dcc, html

from PyMieSim.gui.schemas import FieldSpec, SECTION_FIELDS


def create_layout(default_measure_options: list[str]):
    """Create the top-level Dash layout."""
    return html.Div(
        className="app-shell",
        children=[
            dcc.Store(id="experiment-result"),
            dcc.Download(id="csv-download"),
            html.Header(
                className="hero",
                children=[
                    html.Div(
                        className="hero-copy",
                        children=[
                            html.P("PyMieSim", className="eyebrow"),
                            html.H1("Experiment Lab"),
                            html.P(
                                "Configure source, scatterer, and detector sets, then run the compiled experiment engine directly from Dash.",
                                className="hero-text",
                            ),
                        ],
                    ),
                    html.Div(
                        className="hero-meta",
                        children=[
                            html.Div(className="meta-chip", children=[html.Span("Set-driven"), html.Small("GaussianSet, SphereSet, PhotodiodeSet and more")]),
                            html.Div(className="meta-chip", children=[html.Span("Compiled backend"), html.Small("Sends your configured sets straight into the experiment engine and returns structured results for plotting and export.")]),
                        ],
                    ),
                ],
            ),
            html.Main(
                className="workspace",
                children=[
                    html.Section(
                        className="control-column",
                        children=[
                            html.Section(
                                className="panel plot-controls-panel",
                                children=[
                                    html.Div(className="panel-header", children=[html.H2("Plot Controls")]),
                                    html.P(
                                        "Choose the response to compute and the parameter to place on the horizontal axis.",
                                        className="helper-copy",
                                    ),
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
                                                    ),
                                                ],
                                            ),
                                            html.Div(
                                                className="field-block",
                                                children=[
                                                    html.Label("X Axis", htmlFor="x-axis-select"),
                                                    dcc.Dropdown(id="x-axis-select", className="dashboard-dropdown", options=[], placeholder="Detected from fields with multiple values"),
                                                ],
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                            html.Section(
                                className="panel run-panel",
                                children=[
                                    html.Div(className="panel-header run-panel-header", children=[html.H2("Run Experiment")]),
                                    html.P(
                                        "Launch the compiled experiment engine with the current sets and refresh the figure instantly.",
                                        className="run-panel-copy",
                                    ),
                                    html.Div(
                                        className="run-actions",
                                        children=[
                                            html.Button("Run Experiment", id="run-experiment", n_clicks=0, className="run-button run-button-primary"),
                                            html.Button("Export CSV", id="export-csv", n_clicks=0, className="run-button export-button"),
                                        ],
                                    ),
                                    html.Div(id="status-banner", className="status-banner idle", children="Ready."),
                                ],
                            ),
                            html.Section(
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
                    html.Section(
                        className="result-column",
                        children=[
                            html.Div(
                                className="set-panel-grid",
                                children=[
                                    _section_shell(
                                        "Source Set",
                                        dcc.Dropdown(
                                            id="source-type",
                                            className="dashboard-dropdown",
                                            options=[{"label": key.replace("Set", ""), "value": key} for key in SECTION_FIELDS["source"]],
                                            value="GaussianSet",
                                            clearable=False,
                                        ),
                                        html.Div(id="source-fields"),
                                    ),
                                    _section_shell(
                                        "Scatterer Set",
                                        dcc.Dropdown(
                                            id="scatterer-type",
                                            className="dashboard-dropdown",
                                            options=[{"label": key.replace("Set", ""), "value": key} for key in SECTION_FIELDS["scatterer"]],
                                            value="SphereSet",
                                            clearable=False,
                                        ),
                                        html.Div(id="scatterer-fields"),
                                    ),
                                    _section_shell(
                                        "Detector Set",
                                        dcc.Dropdown(
                                            id="detector-type",
                                            className="dashboard-dropdown",
                                            options=[{"label": "No detector" if key == "None" else key.replace("Set", ""), "value": key} for key in SECTION_FIELDS["detector"]],
                                            value="PhotodiodeSet",
                                            clearable=False,
                                        ),
                                        html.Div(id="detector-fields"),
                                    ),
                                ],
                            ),
                            html.Div(id="summary-cards", className="summary-grid"),
                            html.Section(
                                className="panel graph-panel",
                                children=[dcc.Graph(id="result-graph", config={"displaylogo": False})],
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


def render_fields(section: str, section_type: str):
    """Render all fields for one selected section type."""
    field_specs = SECTION_FIELDS[section][section_type]

    if not field_specs:
        return html.Div("This experiment runs without a detector set.", className="empty-section")

    return html.Div(className="field-grid", children=[render_field(section, field_spec) for field_spec in field_specs])


def render_field(section: str, field_spec: FieldSpec):
    """Render one schema field as a labeled text input."""
    description = field_spec.help_text or _build_default_help_text(field_spec)

    return html.Div(
        className="field-block",
        children=[
            html.Label(field_spec.label),
            dcc.Input(
                id={"kind": "field", "section": section, "name": field_spec.name},
                type="text",
                value=field_spec.default,
                debounce=True,
                placeholder=field_spec.placeholder,
                className="field-input",
            ),
            html.Small(description, className="field-help"),
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