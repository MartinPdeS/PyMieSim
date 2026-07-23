"""Dash interface for the PyMieSim experiment dashboard."""

from __future__ import annotations

import logging
from pathlib import Path
import webbrowser
from typing import Any

from dash import ALL, MATCH, Dash, Input, Output, State, dcc, no_update

from PyMieSim.gui.layout import THEME_DARK, THEME_LIGHT, build_page_with_footer, create_layout, render_fields
from PyMieSim.gui.defaults import DEFAULT_APPLICATION_SETTINGS, DEFAULT_PARTICLE_PLOT_SETTINGS, DEFAULT_SWEEP_PLOT_SETTINGS
from PyMieSim.gui.pages.documentation import build_documentation_page
from PyMieSim.gui.pages.citation import build_citation_page
from PyMieSim.gui.pages.home import build_home_page
from PyMieSim.gui.pages.install_local import build_install_local_page
from PyMieSim.gui.pages.experiment import build_experiment_page
from PyMieSim.gui.pages.settings import build_settings_page
from PyMieSim.gui.pages.single import build_single_page
from PyMieSim.gui.schemas import DETECTOR_FIELDS, SCATTERER_FIELDS, SINGLE_SCATTERER_FIELDS, SINGLE_SOURCE_FIELDS, SOURCE_FIELDS
from PyMieSim.gui.services import (
    available_measures,
    apply_plot_settings,
    build_figure,
    export_result_to_csv,
    export_single_result_to_csv,
    infer_variable_fields,
    build_single_figure,
    run_experiment,
    _parse_field_value,
    _validate_positive_field,
    _is_angular_unit,
)


LOGGER = logging.getLogger(__name__)


def create_dash_app() -> Dash:
    """Create and configure the experiment dashboard Dash application."""
    LOGGER.debug("Creating Dash application")
    app = Dash(
        __name__,
        title="PyMieSim Parameter Sweep Lab",
        assets_folder=str(Path(__file__).with_name("assets")),
        suppress_callback_exceptions=True,
    )
    app.index_string = app.index_string.replace(
        "{%favicon%}",
        '<link rel="icon" type="image/svg+xml" href="/assets/pymiesim-favicon.svg?v=2">'
        '<link rel="shortcut icon" type="image/svg+xml" href="/assets/pymiesim-favicon.svg?v=2">',
    )

    initial_measures = available_measures("SphereSet", "PhotodiodeSet")
    LOGGER.debug("Initial measures loaded: %s", initial_measures)
    app.layout = create_layout(initial_measures)
    _register_callbacks(app, initial_measures)
    LOGGER.debug("Dash application initialized with %d callbacks", len(app.callback_map))
    return app


class OpticalSetupGUI:
    """Backward-compatible wrapper around the new Dash experiment dashboard."""

    def __init__(self):
        """Initialize the Dash application wrapper."""
        LOGGER.debug("Initializing OpticalSetupGUI wrapper")
        self.app = app

    def run(
        self,
        host: str = "0.0.0.0",
        port: str = "8050",
        open_browser: bool = False,
        debug: bool = False,
    ):
        """Run the dashboard server."""
        url = f"http://{host}:{port}/"
        LOGGER.debug("Running dashboard at %s debug=%s open_browser=%s", url, debug, open_browser)

        if open_browser:
            webbrowser.open(url, new=2)

        self.app.run(debug=debug, host=host, port=port)


def _register_callbacks(app: Dash, default_measure_options: list[str]) -> None:
    """Register all dashboard callbacks."""
    LOGGER.debug("Registering dashboard callbacks")

    @app.callback(
        Output("theme-link", "href"),
        Output("theme-store", "data"),
        Output("settings-theme-mode", "value"),
        Output("sidebar-logo", "src"),
        Input("settings-theme-mode", "value"),
        State("theme-store", "data"),
    )
    def _sync_theme(settings_theme: str | None, stored_theme: dict | None):
        theme_mode = settings_theme
        if theme_mode not in {"light", "dark"}:
            theme_mode = (stored_theme or {}).get("theme", DEFAULT_APPLICATION_SETTINGS["theme"])
        mode = "light" if theme_mode == "light" else "dark"
        logo = "/assets/pymiesim-logo.svg" if mode == "light" else "/assets/pymiesim-logo-dark.svg"
        return (THEME_LIGHT if mode == "light" else THEME_DARK), {"theme": mode}, mode, logo

    @app.callback(
        Output("page-content", "children"),
        Output("sidebar-link-home", "className"),
        Output("sidebar-link-experiment", "className"),
        Output("sidebar-link-single", "className"),
        Output("sidebar-link-documentation", "className"),
        Output("sidebar-link-settings", "className"),
        Output("home-visit-count", "data"),
        Input("url", "pathname"),
        State("home-visit-count", "data"),
        State("experiment-run-count", "data"),
        State("single-run-count", "data"),
        State("theme-store", "data"),
        State("plot-settings-store", "data"),
    )
    def _route_pages(pathname: str | None, home_visits: int, experiment_runs: int, single_runs: int, theme_store: dict | None, plot_settings: dict | None):
        """Render only the selected route page inside the persistent shell."""
        route = pathname or "/"
        active = {
            "home": "sidebar-link",
            "experiment": "sidebar-link",
            "single": "sidebar-link",
            "documentation": "sidebar-link",
            "settings": "sidebar-link",
        }
        home_visits = int(home_visits or 0)
        if route == "/":
            home_visits += 1
        metrics = {
            "home_page_visits": home_visits,
            "experiment_runs": int(experiment_runs or 0),
            "single_runs": int(single_runs or 0),
        }
        if route == "/documentation":
            active["documentation"] += " active"
            return build_page_with_footer(build_documentation_page()), *(active[key] for key in ("home", "experiment", "single", "documentation", "settings")), home_visits
        if route == "/citation":
            active["home"] += " active"
            return build_page_with_footer(build_citation_page()), *(active[key] for key in ("home", "experiment", "single", "documentation", "settings")), home_visits
        if route == "/documentation/install-local":
            active["documentation"] += " active"
            return build_page_with_footer(build_install_local_page()), *(active[key] for key in ("home", "experiment", "single", "documentation", "settings")), home_visits
        if route == "/settings":
            active["settings"] += " active"
            return build_page_with_footer(build_settings_page((theme_store or {}).get("theme", "light"), plot_settings or {})), *(active[key] for key in ("home", "experiment", "single", "documentation", "settings")), home_visits
        if route == "/single":
            active["single"] += " active"
            return build_page_with_footer(build_single_page((plot_settings or {}).get("particle_explorer", {}))), *(active[key] for key in ("home", "experiment", "single", "documentation", "settings")), home_visits
        if route == "/experiment":
            active["experiment"] += " active"
            return build_page_with_footer(build_experiment_page(default_measure_options, plot_settings or {})), *(active[key] for key in ("home", "experiment", "single", "documentation", "settings")), home_visits
        active["home"] += " active"
        return build_page_with_footer(build_home_page(metrics)), *(active[key] for key in ("home", "experiment", "single", "documentation", "settings")), home_visits

    @app.callback(
        Output("plot-settings-store", "data"),
        Input("settings-particle-font-size", "value", allow_optional=True),
        Input("settings-particle-line-width", "value", allow_optional=True),
        Input("settings-particle-graph-height", "value", allow_optional=True),
        Input("settings-particle-template", "value", allow_optional=True),
        Input("settings-particle-coordinates", "value", allow_optional=True),
        Input("settings-particle-x-scale", "value", allow_optional=True),
        Input("settings-particle-legend", "value", allow_optional=True),
        Input("settings-particle-grid", "value", allow_optional=True),
        Input("settings-particle-log-y", "value", allow_optional=True),
        Input("settings-sweep-font-size", "value", allow_optional=True),
        Input("settings-sweep-line-width", "value", allow_optional=True),
        Input("settings-sweep-graph-height", "value", allow_optional=True),
        Input("settings-sweep-template", "value", allow_optional=True),
        Input("settings-sweep-coordinates", "value", allow_optional=True),
        Input("settings-sweep-x-scale", "value", allow_optional=True),
        Input("settings-sweep-legend", "value", allow_optional=True),
        Input("settings-sweep-grid", "value", allow_optional=True),
        Input("settings-sweep-log-y", "value", allow_optional=True),
        Input("plot-single-x-scale", "value", allow_optional=True),
        Input("plot-single-y-scale", "value", allow_optional=True),
        Input("plot-single-font-size", "value", allow_optional=True),
        Input("plot-single-line-width", "value", allow_optional=True),
        Input("plot-single-legend", "value", allow_optional=True),
        Input("plot-single-grid", "value", allow_optional=True),
        Input("plot-experiment-x-scale", "value", allow_optional=True),
        Input("plot-experiment-y-scale", "value", allow_optional=True),
        Input("plot-experiment-font-size", "value", allow_optional=True),
        Input("plot-experiment-line-width", "value", allow_optional=True),
        Input("plot-experiment-legend", "value", allow_optional=True),
        Input("plot-experiment-grid", "value", allow_optional=True),
        State("plot-settings-store", "data"),
    )
    def _sync_plot_settings(*values):
        stored_settings = values[-1] or {}
        particle_values = values[:9]
        sweep_values = values[9:18]
        single_local_values = values[18:24]
        experiment_local_values = values[24:30]

        def merge(defaults, key, settings_values):
            current = {**defaults, **(stored_settings.get(key, {}) if isinstance(stored_settings, dict) else {})}
            names = ("font_size", "line_width", "graph_height", "template", "coordinate_system", "x_scale", "show_legend", "show_grid", "log_y")
            for name, value in zip(names, settings_values):
                if value is not None:
                    current[name] = value
            return current

        def merge_local(current, local_values):
            names = ("x_scale", "log_y", "font_size", "line_width", "show_legend", "show_grid")
            for name, value in zip(names, local_values):
                if value is not None:
                    current[name] = value
            return current

        particle = merge(DEFAULT_PARTICLE_PLOT_SETTINGS, "particle_explorer", particle_values)
        sweep = merge(DEFAULT_SWEEP_PLOT_SETTINGS, "parameter_sweep", sweep_values)
        return {
            "particle_explorer": merge_local(particle, single_local_values),
            "parameter_sweep": merge_local(sweep, experiment_local_values),
        }

    @app.callback(
        Output({"kind": "field", "section": MATCH, "name": ALL}, "className"),
        Input({"kind": "field", "section": MATCH, "name": ALL}, "value"),
        State({"kind": "field", "section": MATCH, "name": ALL}, "id"),
        State("source-type", "value", allow_optional=True),
        State("scatterer-type", "value", allow_optional=True),
        State("detector-type", "value", allow_optional=True),
        State("single-source-type", "value", allow_optional=True),
        State("single-scatterer-type", "value", allow_optional=True),
    )
    def _validate_fields(values, field_ids, source_type, scatterer_type, detector_type, single_source_type, single_scatterer_type):
        """Mark schema-invalid dynamic inputs without interrupting the form."""
        if not field_ids:
            return []

        selected_types = {
            "source": source_type,
            "scatterer": scatterer_type,
            "detector": detector_type,
            "single-source": single_source_type,
            "single-scatterer": single_scatterer_type,
        }
        schema_groups = {
            "source": SOURCE_FIELDS,
            "scatterer": SCATTERER_FIELDS,
            "detector": DETECTOR_FIELDS,
            "single-source": SINGLE_SOURCE_FIELDS,
            "single-scatterer": SINGLE_SCATTERER_FIELDS,
        }

        classes = []
        for field_id, raw_value in zip(field_ids, values):
            section = field_id["section"]
            selected_type = selected_types.get(section)
            field_specs = schema_groups.get(section, {}).get(selected_type, ())
            field = next((spec for spec in field_specs if spec.name == field_id["name"]), None)
            valid = field is not None

            if field is not None:
                empty = raw_value is None or str(raw_value).strip() == ""
                valid = field.optional and empty
                if not empty:
                    try:
                        parsed_value = _parse_field_value(field.kind, raw_value, field.unit)
                        if field.name in {"wavelength", "optical_power", "numerical_aperture", "amplitude", "diameter", "core_diameter", "shell_thickness", "sampling"}:
                            _validate_positive_field(field.name, parsed_value)
                        valid = True
                    except (TypeError, ValueError, KeyError):
                        valid = False

            is_detector_default = section == "detector" and field is not None and field.name in {"polarization_filter", "medium"} and (raw_value is None or str(raw_value).strip() in {"", field.default})
            classes.append("field-input field-input-default" if valid and is_detector_default else "field-input" if valid else "field-input field-input-invalid")

        return classes

    @app.callback(Output("source-fields", "children"), Input("source-type", "value"))
    def _render_source_fields(source_type: str):
        LOGGER.debug("Rendering source fields for %s", source_type)
        return render_fields("source", source_type)

    @app.callback(Output("scatterer-fields", "children"), Input("scatterer-type", "value"))
    def _render_scatterer_fields(scatterer_type: str):
        LOGGER.debug("Rendering scatterer fields for %s", scatterer_type)
        return render_fields("scatterer", scatterer_type)

    @app.callback(Output("detector-fields", "children"), Input("detector-type", "value"))
    def _render_detector_fields(detector_type: str):
        LOGGER.debug("Rendering detector fields for %s", detector_type)
        return render_fields("detector", detector_type)

    @app.callback(Output("single-source-fields", "children"), Input("single-source-type", "value"))
    def _render_single_source_fields(source_type: str):
        return render_fields("single-source", source_type)

    @app.callback(Output("single-scatterer-fields", "children"), Input("single-scatterer-type", "value"))
    def _render_single_scatterer_fields(scatterer_type: str):
        return render_fields("single-scatterer", scatterer_type)

    @app.callback(
        Output("measure-select", "options"),
        Output("measure-select", "value"),
        Input("scatterer-type", "value"),
        Input("detector-type", "value"),
        State("measure-select", "value"),
    )
    def _update_measure_options(scatterer_type: str, detector_type: str, current_measure: str | None):
        LOGGER.debug(
            "Updating measure options for scatterer=%s detector=%s current=%s",
            scatterer_type,
            detector_type,
            current_measure,
        )
        measures = available_measures(scatterer_type, detector_type)
        options = [{"label": measure, "value": measure} for measure in measures]
        value = current_measure if current_measure in measures else measures[0]
        return options, value

    @app.callback(
        Output("x-axis-select", "options"),
        Output("x-axis-select", "value"),
        Input("source-type", "value"),
        Input({"kind": "field", "section": "source", "name": ALL}, "value"),
        Input({"kind": "field", "section": "source", "name": ALL}, "id"),
        Input("scatterer-type", "value"),
        Input({"kind": "field", "section": "scatterer", "name": ALL}, "value"),
        Input({"kind": "field", "section": "scatterer", "name": ALL}, "id"),
        Input("detector-type", "value"),
        Input({"kind": "field", "section": "detector", "name": ALL}, "value"),
        Input({"kind": "field", "section": "detector", "name": ALL}, "id"),
        State("x-axis-select", "value"),
    )
    def _update_x_axis_options(
        source_type: str,
        source_values: list[str],
        source_ids: list[dict[str, str]],
        scatterer_type: str,
        scatterer_values: list[str],
        scatterer_ids: list[dict[str, str]],
        detector_type: str,
        detector_values: list[str],
        detector_ids: list[dict[str, str]],
        current_x_axis: str | None,
    ):
        LOGGER.debug(
            "Updating x-axis options from form state source=%s scatterer=%s detector=%s current=%s",
            source_type,
            scatterer_type,
            detector_type,
            current_x_axis,
        )

        variable_fields = infer_variable_fields(
            source_type=source_type,
            source_values=_pair_ids_with_values(source_ids, source_values),
            scatterer_type=scatterer_type,
            scatterer_values=_pair_ids_with_values(scatterer_ids, scatterer_values),
            detector_type=detector_type,
            detector_values=_pair_ids_with_values(detector_ids, detector_values),
        )

        options = [{"label": field_name, "value": field_name} for field_name in variable_fields]
        value = current_x_axis if current_x_axis in variable_fields else (variable_fields[0] if variable_fields else None)
        return options, value

    @app.callback(
        Output("experiment-result", "data"),
        Output("experiment-run-count", "data"),
        Input("source-type", "value"),
        Input({"kind": "field", "section": "source", "name": ALL}, "value"),
        State({"kind": "field", "section": "source", "name": ALL}, "id"),
        Input("scatterer-type", "value"),
        Input({"kind": "field", "section": "scatterer", "name": ALL}, "value"),
        State({"kind": "field", "section": "scatterer", "name": ALL}, "id"),
        Input("detector-type", "value"),
        Input({"kind": "field", "section": "detector", "name": ALL}, "value"),
        State({"kind": "field", "section": "detector", "name": ALL}, "id"),
        Input("measure-select", "value"),
        State("experiment-run-count", "data"),
        prevent_initial_call=True,
    )
    def _run_experiment(
        source_type: str,
        source_values: list[str],
        source_ids: list[dict[str, str]],
        scatterer_type: str,
        scatterer_values: list[str],
        scatterer_ids: list[dict[str, str]],
        detector_type: str,
        detector_values: list[str],
        detector_ids: list[dict[str, str]],
        measure: str,
        experiment_runs: int,
    ):
        next_experiment_runs = int(experiment_runs or 0)

        LOGGER.debug(
            "Auto-running parameter sweep with source=%s scatterer=%s detector=%s measure=%s",
            source_type,
            scatterer_type,
            detector_type,
            measure,
        )

        try:
            result = run_experiment(
                source_type=source_type,
                source_values=_pair_ids_with_values(source_ids, source_values),
                scatterer_type=scatterer_type,
                scatterer_values=_pair_ids_with_values(scatterer_ids, scatterer_values),
                detector_type=detector_type,
                detector_values=_pair_ids_with_values(detector_ids, detector_values),
                measure=measure,
            )
        except Exception as error:
            LOGGER.debug("Skipping incomplete or invalid auto-run input: %s", error)
            return no_update, next_experiment_runs

        LOGGER.debug(
            "Experiment run finished rows=%d x_axis_options=%s",
            result["row_count"],
            result["parameter_columns"],
        )

        return result, next_experiment_runs + 1

    @app.callback(
        Output("csv-download", "data"),
        Input("export-csv", "n_clicks"),
        State("experiment-result", "data"),
        State("measure-select", "value"),
        prevent_initial_call=True,
    )
    def _export_csv(n_clicks: int, result: dict | None, measure: str | None):
        LOGGER.debug("CSV export requested for measure=%s", measure)

        if not n_clicks:
            LOGGER.debug("Ignoring CSV callback without an explicit button click")
            return no_update

        csv_content = export_result_to_csv(result)

        if not csv_content:
            LOGGER.debug("No CSV exported because result content is empty")
            return no_update

        filename = f"pymiesim_{measure or 'result'}.csv"
        LOGGER.debug("CSV export generated filename=%s", filename)
        return dcc.send_string(csv_content, filename)

    @app.callback(
        Output("single-csv-download", "data"),
        Input("export-single-csv", "n_clicks"),
        State("single-result", "data"),
        prevent_initial_call=True,
    )
    def _export_single_csv(n_clicks: int, result: dict | None):
        if not n_clicks:
            return no_update

        csv_content = export_single_result_to_csv(result)
        if not csv_content:
            return no_update

        return dcc.send_string(csv_content, "pymiesim_particle_explorer.csv")

    @app.callback(
        Output("result-graph", "figure"),
        Input("experiment-result", "data"),
        Input("x-axis-select", "value"),
        Input("plot-settings-store", "data"),
        Input("theme-store", "data"),
        Input("plot-experiment-x-scale", "value", allow_optional=True),
        Input("plot-experiment-y-scale", "value", allow_optional=True),
        Input("plot-experiment-font-size", "value", allow_optional=True),
        Input("plot-experiment-line-width", "value", allow_optional=True),
        Input("plot-experiment-legend", "value", allow_optional=True),
        Input("plot-experiment-grid", "value", allow_optional=True),
        Input("plot-experiment-projection", "value", allow_optional=True),
    )
    def _update_outputs(result: dict | None, x_axis: str | None, plot_settings: dict | None, theme_store: dict | None, *local_values):
        LOGGER.debug("Updating outputs for x_axis=%s result_present=%s", x_axis, bool(result))
        projection = local_values[-1] or "cartesian"
        sweep_settings = (plot_settings or {}).get("parameter_sweep", plot_settings or {})
        sweep_settings = _merge_local_plot_values(sweep_settings, local_values[:-1])
        return build_figure(result, x_axis, plot_settings=sweep_settings, theme=(theme_store or {}).get("theme", "light"), projection=projection)

    @app.callback(
        Output("single-result", "data"),
        Output("single-run-count", "data"),
        Input("single-source-type", "value"),
        Input({"kind": "field", "section": "single-source", "name": ALL}, "value"),
        State({"kind": "field", "section": "single-source", "name": ALL}, "id"),
        Input("single-scatterer-type", "value"),
        Input({"kind": "field", "section": "single-scatterer", "name": ALL}, "value"),
        State({"kind": "field", "section": "single-scatterer", "name": ALL}, "id"),
        Input("single-representation", "value"),
        Input("single-projection", "value"),
        Input("single-sampling", "value"),
        State("single-run-count", "data"),
    )
    def _run_single(
        source_type: str,
        source_values: list[str],
        source_ids: list[dict[str, str]],
        scatterer_type: str,
        scatterer_values: list[str],
        scatterer_ids: list[dict[str, str]],
        representation: str,
        projection: str,
        sampling: int,
        single_runs: int,
    ):
        next_single_runs = int(single_runs or 0)
        try:
            figure, summary = build_single_figure(
                source_type=source_type,
                source_values=_pair_ids_with_values(source_ids, source_values),
                scatterer_type=scatterer_type,
                scatterer_values=_pair_ids_with_values(scatterer_ids, scatterer_values),
                representation=representation,
                projection=projection,
                sampling=sampling or 120,
            )
        except Exception as error:
            LOGGER.exception("Single representation render failed")
            return None, next_single_runs

        return {"figure": figure.to_plotly_json(), "summary": summary}, next_single_runs + 1

    @app.callback(
        Output("single-projection", "options"),
        Output("single-projection", "value"),
        Input("single-representation", "value"),
        State("single-projection", "value"),
    )
    def _update_single_projection_options(representation: str, current_projection: str | None):
        if representation == "s1s2":
            options = [{"label": "2D plot", "value": "2d"}, {"label": "Polar 1D", "value": "polar_1d"}]
            return options, current_projection if current_projection == "polar_1d" else "2d"

        supports_3d = representation in {"stokes", "stokes_q", "stokes_u", "stokes_v", "spf", "farfields"}
        options = [{"label": "2D heatmap" if supports_3d else "2D plot", "value": "2d"}]
        if supports_3d:
            options.extend([
                {"label": "3D sphere", "value": "3d"},
                {"label": "3D radial surface", "value": "3d_radial"},
            ])
        return options, current_projection if supports_3d and current_projection in {"3d", "3d_radial"} else "2d"

    @app.callback(
        Output("single-plot-options-container", "children"),
        Input("single-representation", "value"),
        Input("single-projection", "value"),
        State("plot-settings-store", "data"),
    )
    def _update_single_plot_options(representation: str, projection: str, plot_settings: dict | None):
        particle_settings = (plot_settings or {}).get("particle_explorer", plot_settings or {})
        from PyMieSim.gui.layout import _plot_options_card

        return _plot_options_card("single", particle_settings, representation, projection)

    @app.callback(
        Output("experiment-plot-options-container", "children"),
        Input("x-axis-select", "value"),
        Input("experiment-result", "data"),
        State("plot-settings-store", "data"),
        State("plot-experiment-projection", "value", allow_optional=True),
    )
    def _update_experiment_plot_options(x_axis: str | None, result: dict | None, plot_settings: dict | None, current_projection: str | None):
        sweep_settings = (plot_settings or {}).get("parameter_sweep", plot_settings or {})
        units = (result or {}).get("units", {})
        xaxis_unit = units.get(x_axis)
        if xaxis_unit is None and x_axis:
            matching_units = [value for key, value in units.items() if key.rsplit(":", 1)[-1] == x_axis]
            xaxis_unit = matching_units[0] if matching_units else None
        is_angle = _is_angular_unit(xaxis_unit)
        options = [{"label": "Cartesian", "value": "cartesian"}]
        if is_angle:
            options.append({"label": "Polar", "value": "polar"})
        projection = current_projection if is_angle and current_projection == "polar" else "cartesian"
        from PyMieSim.gui.layout import _plot_options_card

        return _plot_options_card("experiment", sweep_settings, projection=projection, projection_options=options)

    @app.callback(
        Output("single-graph", "figure"),
        Input("single-result", "data"),
        Input("plot-settings-store", "data"),
        Input("theme-store", "data"),
        Input("plot-single-x-scale", "value", allow_optional=True),
        Input("plot-single-y-scale", "value", allow_optional=True),
        Input("plot-single-font-size", "value", allow_optional=True),
        Input("plot-single-line-width", "value", allow_optional=True),
        Input("plot-single-legend", "value", allow_optional=True),
        Input("plot-single-grid", "value", allow_optional=True),
    )
    def _update_single_outputs(result: dict | None, plot_settings: dict | None, theme_store: dict | None, *local_values):
        particle_settings = (plot_settings or {}).get("particle_explorer", plot_settings or {})
        particle_settings = _merge_local_plot_values(particle_settings, local_values)
        if not result:
            return build_single_empty_figure(plot_settings=particle_settings, theme=(theme_store or {}).get("theme", "light"))
        from plotly.graph_objects import Figure

        return apply_plot_settings(Figure(result["figure"]), particle_settings, (theme_store or {}).get("theme", "light"))


def _pair_ids_with_values(ids: list[dict[str, str]], values: list[str]) -> dict[str, str]:
    """Convert dynamic Dash field IDs and values into a flat mapping."""
    return {field_id["name"]: value for field_id, value in zip(ids, values)}


def _merge_local_plot_values(settings: dict | None, values: tuple[Any, ...]) -> dict:
    """Overlay values from a plot-local options card onto stored preferences."""
    merged = dict(settings or {})
    names = ("x_scale", "log_y", "font_size", "line_width", "show_legend", "show_grid")
    for name, value in zip(names, values):
        if value is not None:
            merged[name] = value
    return merged


def build_single_empty_figure(plot_settings: dict | None = None, theme: str = "light"):
    """Return the initial representation canvas."""
    from plotly.graph_objects import Figure

    figure = Figure()
    figure.update_layout(template="plotly_white", xaxis_visible=False, yaxis_visible=False, annotations=[{"text": "Configure a setup to inspect a representation.", "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5, "showarrow": False}])
    return apply_plot_settings(figure, plot_settings, theme)


app = create_dash_app()
server = app.server


if __name__ == "__main__":
    _gui = OpticalSetupGUI()
    _gui.run(host="0.0.0.0", port="8050", open_browser=False, debug=True)
