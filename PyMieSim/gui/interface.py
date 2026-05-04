"""Dash interface for the PyMieSim experiment dashboard."""

from __future__ import annotations

import logging
from pathlib import Path
import webbrowser

from dash import ALL, Dash, Input, Output, State, dcc, no_update

from PyMieSim.gui.layout import create_layout, render_fields, render_summary_cards
from PyMieSim.gui.services import (
    available_measures,
    build_figure,
    build_summary,
    export_result_to_csv,
    infer_variable_fields,
    run_experiment,
)


LOGGER = logging.getLogger(__name__)


def create_dash_app() -> Dash:
    """Create and configure the experiment dashboard Dash application."""
    LOGGER.debug("Creating Dash application")
    app = Dash(
        __name__,
        title="PyMieSim Experiment Lab",
        assets_folder=str(Path(__file__).with_name("assets")),
        suppress_callback_exceptions=True,
    )

    initial_measures = available_measures("SphereSet", "PhotodiodeSet")
    LOGGER.debug("Initial measures loaded: %s", initial_measures)
    app.layout = create_layout(initial_measures)
    _register_callbacks(app)
    LOGGER.debug("Dash application initialized with %d callbacks", len(app.callback_map))
    return app


class OpticalSetupGUI:
    """Backward-compatible wrapper around the new Dash experiment dashboard."""

    def __init__(self):
        """Initialize the Dash application wrapper."""
        LOGGER.debug("Initializing OpticalSetupGUI wrapper")
        self.app = create_dash_app()

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


def _register_callbacks(app: Dash) -> None:
    """Register all dashboard callbacks."""
    LOGGER.debug("Registering dashboard callbacks")

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
        Output("status-banner", "children"),
        Output("status-banner", "className"),
        Input("run-experiment", "n_clicks"),
        State("source-type", "value"),
        State({"kind": "field", "section": "source", "name": ALL}, "value"),
        State({"kind": "field", "section": "source", "name": ALL}, "id"),
        State("scatterer-type", "value"),
        State({"kind": "field", "section": "scatterer", "name": ALL}, "value"),
        State({"kind": "field", "section": "scatterer", "name": ALL}, "id"),
        State("detector-type", "value"),
        State({"kind": "field", "section": "detector", "name": ALL}, "value"),
        State({"kind": "field", "section": "detector", "name": ALL}, "id"),
        State("measure-select", "value"),
        prevent_initial_call=True,
    )
    def _run_experiment(
        n_clicks: int,
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
    ):
        del n_clicks

        LOGGER.debug(
            "Run button pressed with source=%s scatterer=%s detector=%s measure=%s",
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
            LOGGER.exception("Experiment run failed")
            return None, str(error), "status-banner error"

        LOGGER.debug(
            "Experiment run finished rows=%d x_axis_options=%s",
            result["row_count"],
            result["parameter_columns"],
        )

        return result, "Experiment completed. Figure updated and CSV is ready.", "status-banner success"

    @app.callback(
        Output("csv-download", "data"),
        Input("export-csv", "n_clicks"),
        State("experiment-result", "data"),
        State("measure-select", "value"),
        prevent_initial_call=True,
    )
    def _export_csv(n_clicks: int, result: dict | None, measure: str | None):
        del n_clicks
        LOGGER.debug("CSV export requested for measure=%s", measure)

        csv_content = export_result_to_csv(result)

        if not csv_content:
            LOGGER.debug("No CSV exported because result content is empty")
            return no_update

        filename = f"pymiesim_{measure or 'result'}.csv"
        LOGGER.debug("CSV export generated filename=%s", filename)
        return dcc.send_string(csv_content, filename)

    @app.callback(
        Output("result-graph", "figure"),
        Output("summary-cards", "children"),
        Input("experiment-result", "data"),
        Input("x-axis-select", "value"),
    )
    def _update_outputs(result: dict | None, x_axis: str | None):
        LOGGER.debug("Updating outputs for x_axis=%s result_present=%s", x_axis, bool(result))
        figure = build_figure(result, x_axis)
        summary = render_summary_cards(build_summary(result))

        if not result:
            return figure, summary

        return figure, summary


def _pair_ids_with_values(ids: list[dict[str, str]], values: list[str]) -> dict[str, str]:
    """Convert dynamic Dash field IDs and values into a flat mapping."""
    return {field_id["name"]: value for field_id, value in zip(ids, values)}


if __name__ == "__main__":
    _gui = OpticalSetupGUI()
    _gui.run(host="0.0.0.0", port="8050", open_browser=False, debug=True)
