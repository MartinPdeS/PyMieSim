"""Pure backend services for the experiment dashboard."""

from __future__ import annotations

import logging
from typing import Any, Dict

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from PyMieSim.experiment import Setup
from PyMieSim.experiment.detector_set import CoherentModeSet, PhotodiodeSet
from PyMieSim.experiment.scatterer_set import CoreShellSet, InfiniteCylinderSet, SphereSet
from PyMieSim.experiment.source_set import GaussianSet, PlaneWaveSet
from PyMieSim.single import scatterer as single_scatterer
from PyMieSim.single import source as single_source
from PyMieSim.single import Setup as SingleSetup
from PyMieSim.units import ureg

from PyMieSim.gui.parsing import (
    parse_material_values,
    parse_mode_numbers,
    parse_numeric_expression,
    parse_polarization,
    parse_quantity_expression,
    serialize_value,
)
from PyMieSim.gui.defaults import DEFAULT_PLOT_SETTINGS
from PyMieSim.gui.schemas import (
    DETECTOR_FIELDS,
    SCATTERER_FIELDS,
    SINGLE_SCATTERER_FIELDS,
    SINGLE_SOURCE_FIELDS,
    SOURCE_FIELDS,
)


LOGGER = logging.getLogger(__name__)

_POSITIVE_FIELDS = frozenset({
    "wavelength",
    "optical_power",
    "numerical_aperture",
    "amplitude",
    "diameter",
    "core_diameter",
    "shell_thickness",
    "sampling",
})


SOURCE_TYPES = {
    "GaussianSet": GaussianSet,
    "PlaneWaveSet": PlaneWaveSet,
}

SCATTERER_TYPES = {
    "SphereSet": SphereSet,
    "InfiniteCylinderSet": InfiniteCylinderSet,
    "CoreShellSet": CoreShellSet,
}

DETECTOR_TYPES = {
    "None": None,
    "PhotodiodeSet": PhotodiodeSet,
    "CoherentModeSet": CoherentModeSet,
}

SINGLE_SOURCE_TYPES = {"Gaussian": single_source.Gaussian, "PlaneWave": single_source.PlaneWave}
SINGLE_SCATTERER_TYPES = {
    "Sphere": single_scatterer.Sphere,
    "InfiniteCylinder": single_scatterer.InfiniteCylinder,
    "CoreShell": single_scatterer.CoreShell,
}


def available_measures(scatterer_type: str, detector_type: str) -> list[str]:
    """Return measures valid for the selected experiment configuration."""
    LOGGER.debug("Computing available measures for scatterer=%s detector=%s", scatterer_type, detector_type)
    measures = list(SCATTERER_TYPES[scatterer_type].available_measure_list)

    if detector_type == "None":
        measures = [measure for measure in measures if measure != "coupling"]

    return measures


def build_source_set(source_type: str, raw_values: Dict[str, Any]) -> Any:
    """Build the selected source set from raw dashboard values."""
    LOGGER.debug("Building source set %s with raw values: %s", source_type, raw_values)
    return SOURCE_TYPES[source_type](**_parse_section_fields(SOURCE_FIELDS[source_type], raw_values))


def build_scatterer_set(scatterer_type: str, raw_values: Dict[str, Any]) -> Any:
    """Build the selected scatterer set from raw dashboard values."""
    LOGGER.debug("Building scatterer set %s with raw values: %s", scatterer_type, raw_values)
    return SCATTERER_TYPES[scatterer_type](**_parse_section_fields(SCATTERER_FIELDS[scatterer_type], raw_values))


def build_detector_set(detector_type: str, raw_values: Dict[str, Any]) -> Any:
    """Build the selected detector set from raw dashboard values."""
    LOGGER.debug("Building detector set %s with raw values: %s", detector_type, raw_values)
    detector_class = DETECTOR_TYPES[detector_type]

    if detector_class is None:
        LOGGER.debug("Detector type is None; running detector-free experiment")
        return None

    return detector_class(**_parse_section_fields(DETECTOR_FIELDS[detector_type], raw_values))


def build_single_setup(
    *, source_type: str, source_values: Dict[str, Any], scatterer_type: str, scatterer_values: Dict[str, Any]
) -> Any:
    """Build a single-scatterer setup from the representation tab fields."""
    source = SINGLE_SOURCE_TYPES[source_type](**_parse_section_fields(SINGLE_SOURCE_FIELDS[source_type], source_values))
    scatterer = SINGLE_SCATTERER_TYPES[scatterer_type](**_parse_section_fields(SINGLE_SCATTERER_FIELDS[scatterer_type], scatterer_values))
    return SingleSetup(scatterer=scatterer, source=source)


def build_single_figure(
    *,
    source_type: str,
    source_values: Dict[str, Any],
    scatterer_type: str,
    scatterer_values: Dict[str, Any],
    representation: str,
    sampling: int,
    projection: str = "2d",
    nearfield_mode: str = "absolute",
    plot_settings: dict[str, Any] | None = None,
    theme: str = "light",
) -> tuple[go.Figure, dict[str, str]]:
    """Compute one single-scatterer representation and turn it into Plotly."""
    setup = build_single_setup(
        source_type=source_type,
        source_values=source_values,
        scatterer_type=scatterer_type,
        scatterer_values=scatterer_values,
    )
    sampling = max(24, min(int(sampling), 300))
    figure = go.Figure()

    if representation == "s1s2":
        angles = np.linspace(-np.pi, np.pi, sampling) * ureg.radian
        s1, s2 = setup.get_s1s2(angles=angles)
        x = angles.to(ureg.degree).magnitude
        if projection == "polar_1d":
            figure.add_trace(go.Scatterpolar(theta=x, r=np.abs(s1.magnitude), mode="lines", name="|S1|"))
            figure.add_trace(go.Scatterpolar(theta=x, r=np.abs(s2.magnitude), mode="lines", name="|S2|"))
            figure.update_layout(polar={"angularaxis": {"direction": "counterclockwise", "rotation": 90}, "radialaxis": {"title": {"text": ""}}})
        else:
            figure.add_trace(go.Scatter(x=x, y=np.abs(s1.magnitude), mode="lines", name="|S1|"))
            figure.add_trace(go.Scatter(x=x, y=np.abs(s2.magnitude), mode="lines", name="|S2|"))
            figure.update_xaxes(title="Scattering angle [degree]")
            figure.update_yaxes(title=f"Amplitude [{_unit_label(s1, 'meter')}]")
        figure.update_layout(meta={"polar_axis_unit": "degree"})
        title = "S1 / S2 scattering amplitudes"
    elif representation.startswith("nearfields"):
        nearfield_components = {
            "nearfields": "|E|",
            "nearfields_ex": "Ex",
            "nearfields_ey": "Ey",
            "nearfields_ez": "Ez",
        }
        component = nearfield_components.get(representation, "|E|")
        nearfields = setup.get_representation("nearfields")
        field = nearfields.compute(component, type="total", sampling=sampling)[component]
        complex_values = np.asarray(field, dtype=complex)
        values = np.abs(complex_values) if nearfield_mode != "real" else complex_values.real
        x = np.asarray(nearfields.u.to(ureg.nanometer).magnitude, dtype=float)
        y = np.asarray(nearfields.v.to(ureg.nanometer).magnitude, dtype=float)
        figure = go.Figure(
            data=go.Heatmap(
                x=x,
                y=y,
                z=values,
                colorscale="Viridis",
                colorbar={"title": "V/m"},
            )
        )
        figure.update_layout(yaxis={"scaleanchor": "x", "scaleratio": 1})
        figure.update_xaxes(title="Plane u [nanometer]")
        figure.update_yaxes(title="Plane v [nanometer]")
        value_label = "|" if nearfield_mode != "real" else "Re("
        value_suffix = "|" if nearfield_mode != "real" else ")"
        title = f"Near-field {value_label}{component.strip('|')}{value_suffix}"
    else:
        backend_representation = "stokes" if representation == "stokes" or representation.startswith("stokes_") else representation
        representation_object = setup.get_representation(backend_representation, sampling=sampling)
        if representation == "stokes" or representation.startswith("stokes_"):
            field = representation.rsplit("_", 1)[-1].upper() if "_" in representation else "I"
            quantity = getattr(representation_object, field)
            values = quantity.magnitude
            value_label = f"Stokes {field}"
        elif representation == "spf":
            field = "SPF"
            quantity = representation_object.SPF
            values = getattr(quantity, "magnitude", quantity)
            value_label = "Scattering intensity"
        elif representation == "farfields":
            field = "total_intensity"
            quantity = np.abs(representation_object.E_phi) ** 2 + np.abs(representation_object.E_theta) ** 2
            values = quantity.magnitude
            value_label = "Intensity"
        else:
            raise ValueError(f"Representation '{representation}' is not supported by the dashboard.")

        values = np.asarray(values)
        value_unit = _unit_label(quantity, "dimensionless")
        if values.ndim == 2:
            azimuth = np.linspace(-180.0, 180.0, values.shape[1])
            polar = np.linspace(-90.0, 90.0, values.shape[0])
            if projection in {"3d", "3d_radial"}:
                azimuth_grid, polar_grid = np.meshgrid(np.deg2rad(azimuth), np.deg2rad(polar))
                sphere_x = np.cos(polar_grid) * np.cos(azimuth_grid)
                sphere_y = np.cos(polar_grid) * np.sin(azimuth_grid)
                sphere_z = np.sin(polar_grid)
                surface_values = np.abs(values)
                if projection == "3d_radial":
                    maximum = float(np.nanmax(surface_values)) if surface_values.size else 0.0
                    radius = surface_values / maximum if maximum > 0.0 else np.ones_like(surface_values)
                    surface_x = sphere_x * radius
                    surface_y = sphere_y * radius
                    surface_z = sphere_z * radius
                else:
                    surface_x, surface_y, surface_z = sphere_x, sphere_y, sphere_z
                figure = go.Figure(data=go.Surface(x=surface_x, y=surface_y, z=surface_z, surfacecolor=surface_values, colorscale="Viridis", colorbar={"title": value_unit}))
                figure.update_layout(scene={"xaxis_title": "x", "yaxis_title": "y", "zaxis_title": "z", "aspectmode": "cube"})
            else:
                figure = go.Figure(data=go.Heatmap(x=azimuth, y=polar, z=values, colorscale="Viridis", colorbar={"title": value_unit}))
                figure.update_xaxes(title="Azimuth angle [degree]")
                figure.update_yaxes(title="Polar angle [degree]")
        else:
            figure.add_trace(go.Scatter(y=values.ravel(), mode="lines", name=field))
            figure.update_xaxes(title="Scattering angle [degree]")
            figure.update_yaxes(title=f"{value_label} [{value_unit}]")
        title = f"{value_label} representation" if representation.startswith("stokes") else f"{representation.title()} representation"

    figure.update_layout(
        template="plotly_white",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        plot_bgcolor="white",
        margin={"l": 45, "r": 20, "t": 55, "b": 45},
        title=title,
        legend_title_text="",
    )
    apply_plot_settings(figure, plot_settings, theme)
    return figure, {"Representation": representation.upper(), "Sampling": str(sampling), "Scatterer": scatterer_type}


def _unit_label(value: Any, fallback: str) -> str:
    """Return a compact display label for a Pint quantity's units."""
    units = getattr(value, "units", None)
    return str(units) if units is not None else fallback


def infer_variable_fields(
    *,
    source_type: str,
    source_values: Dict[str, Any],
    scatterer_type: str,
    scatterer_values: Dict[str, Any],
    detector_type: str,
    detector_values: Dict[str, Any],
) -> list[str]:
    """Return field names whose current input expands to more than one value."""
    variable_fields: list[str] = []

    variable_fields.extend(_infer_variable_fields_for_section(SOURCE_FIELDS[source_type], source_values))
    variable_fields.extend(_infer_variable_fields_for_section(SCATTERER_FIELDS[scatterer_type], scatterer_values))

    if detector_type != "None":
        variable_fields.extend(_infer_variable_fields_for_section(DETECTOR_FIELDS[detector_type], detector_values))

    LOGGER.debug("Detected variable fields from current inputs: %s", variable_fields)
    return variable_fields


def run_experiment(
    *,
    source_type: str,
    source_values: Dict[str, Any],
    scatterer_type: str,
    scatterer_values: Dict[str, Any],
    detector_type: str,
    detector_values: Dict[str, Any],
    measure: str,
) -> dict[str, Any]:
    """Execute the selected experiment and serialize the result for Dash."""
    LOGGER.debug(
        "Starting experiment run source=%s scatterer=%s detector=%s measure=%s",
        source_type,
        scatterer_type,
        detector_type,
        measure,
    )
    source_set = build_source_set(source_type, source_values)
    scatterer_set = build_scatterer_set(scatterer_type, scatterer_values)
    detector_set = build_detector_set(detector_type, detector_values)

    LOGGER.debug("Source set built: %r", source_set)
    LOGGER.debug("Scatterer set built: %r", scatterer_set)
    LOGGER.debug("Detector set built: %r", detector_set)

    setup = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set,
        detector_set=detector_set,
    )

    LOGGER.debug("Experiment setup instantiated: %r", setup)

    dataframe = setup.get(measure, drop_unique_level=True)
    LOGGER.debug("Experiment returned dataframe with shape %s", getattr(dataframe, "shape", None))
    frame = pd.DataFrame(dataframe).copy()
    units = {key: str(value) for key, value in dataframe.attrs.get("units", {}).items()}

    for column in frame.columns:
        frame[column] = frame[column].map(serialize_value)

    parameter_columns = [column for column in frame.columns if column != measure]

    return {
        "measure": measure,
        "units": units,
        "parameter_columns": parameter_columns,
        "rows": frame.to_dict("records"),
        "row_count": len(frame),
    }


def export_result_to_csv(result: dict[str, Any] | None) -> str:
    """Serialize a stored result payload to CSV text."""
    if not result or not result.get("rows"):
        LOGGER.debug("CSV export requested without result rows")
        return ""

    frame = pd.DataFrame(result["rows"])
    LOGGER.debug("Exporting %d rows to CSV", len(frame))
    return frame.to_csv(index=False)


def export_single_result_to_csv(result: dict[str, Any] | None) -> str:
    """Serialize a Particle Explorer figure payload to CSV text."""
    if not result or not result.get("figure"):
        LOGGER.debug("Particle Explorer CSV export requested without a result")
        return ""

    traces = result["figure"].get("data", [])
    rows: list[dict[str, Any]] = []
    for trace in traces:
        name = trace.get("name", "value")
        if "z" in trace:
            z_values = np.asarray(trace["z"])
            x_values = trace.get("x") or list(range(z_values.shape[1]))
            y_values = trace.get("y") or list(range(z_values.shape[0]))
            for y_index, y_value in enumerate(y_values):
                for x_index, x_value in enumerate(x_values):
                    rows.append({"series": name, "x": x_value, "y": y_value, "value": z_values[y_index, x_index]})
            continue

        y_values = trace.get("y", [])
        x_values = trace.get("x") or list(range(len(y_values)))
        rows.extend({"series": name, "x": x_value, "value": y_value} for x_value, y_value in zip(x_values, y_values))

    return pd.DataFrame(rows).to_csv(index=False) if rows else ""


def build_figure(result: dict[str, Any] | None, x_axis: str | None, plot_settings: dict[str, Any] | None = None, theme: str = "light", projection: str = "cartesian") -> go.Figure:
    """Build a Plotly figure from serialized experiment results."""
    LOGGER.debug("Building figure for x_axis=%s", x_axis)
    if not result or not result.get("rows"):
        LOGGER.debug("No result rows available; returning empty figure")
        figure = go.Figure()
        figure.update_layout(
            template="plotly_white",
            xaxis_visible=False,
            yaxis_visible=False,
            annotations=[
                {
                    "text": "Run an experiment to visualize the result.",
                    "xref": "paper",
                    "yref": "paper",
                    "x": 0.5,
                    "y": 0.5,
                    "showarrow": False,
                    "font": {"size": 16},
                }
            ],
        )
        apply_plot_settings(figure, plot_settings, theme)
        return figure

    frame = pd.DataFrame(result["rows"])
    measure = result["measure"]
    parameter_columns = result["parameter_columns"]
    units = result.get("units", {})
    LOGGER.debug("Figure frame columns=%s row_count=%d", list(frame.columns), len(frame))

    if not parameter_columns:
        resolved_x_axis = "index"
        figure = px.bar(frame.assign(index=[0] * len(frame)), x="index", y=measure)
    else:
        resolved_x_axis = _resolve_x_axis(x_axis, parameter_columns)
        xaxis_unit = units.get(resolved_x_axis)
        if xaxis_unit is None:
            matching_units = [value for key, value in units.items() if key.rsplit(":", 1)[-1] == resolved_x_axis]
            xaxis_unit = matching_units[0] if matching_units else None
        extra_axes = [column for column in parameter_columns if column != resolved_x_axis]
        plot_kwargs: dict[str, Any] = {"x": resolved_x_axis, "y": measure}

        if extra_axes:
            series_column = extra_axes[0]
            frame = frame.copy()

            if len(extra_axes) > 1:
                series_column = "__series__"
                frame[series_column] = frame[extra_axes].apply(
                    lambda row: " | ".join(_format_plot_value(value) for value in row),
                    axis=1,
                )
            else:
                frame[series_column] = frame[series_column].map(
                    _format_plot_value
                )

            plot_kwargs["color"] = series_column
            plot_kwargs["line_group"] = series_column

        if projection == "polar" and _is_angular_unit(xaxis_unit):
            figure = go.Figure()
            groups = frame.groupby(series_column, sort=False) if extra_axes else [(measure, frame)]
            for series_name, group in groups:
                # Plotly's polar theta values are degrees by default.  The
                # experiment data can legitimately be returned in radians,
                # so convert at this boundary instead of letting a 0..2*pi
                # sweep collapse into a tiny line around 0 degrees.
                theta = _polar_theta_degrees(group[resolved_x_axis], xaxis_unit)
                radius = pd.to_numeric(group[measure], errors="coerce").to_numpy(dtype=float)
                valid = np.isfinite(theta) & np.isfinite(radius)
                theta = theta[valid]
                radius = radius[valid]
                order = np.argsort(theta)
                figure.add_trace(go.Scatterpolar(theta=theta[order], r=radius[order], mode="lines", name=str(series_name)))
            figure.update_layout(polar={"angularaxis": {"direction": "counterclockwise", "rotation": 90}})
        else:
            figure = px.line(frame, **plot_kwargs)

    measure_unit = units.get(measure)
    xaxis_unit = locals().get("xaxis_unit", units.get(resolved_x_axis))
    yaxis_title = f"{measure} [{measure_unit}]" if measure_unit else measure
    xaxis_title = f"{resolved_x_axis} [{xaxis_unit}]" if xaxis_unit else resolved_x_axis

    figure.update_layout(
        template="plotly_white",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        plot_bgcolor="white",
        margin={"l": 40, "r": 20, "t": 50, "b": 40},
        title=f"{measure} response",
        xaxis_title=None if projection == "polar" and _is_angular_unit(xaxis_unit) else xaxis_title,
        yaxis_title=None if projection == "polar" and _is_angular_unit(xaxis_unit) else yaxis_title,
        legend_title_text="",
        meta={
            "polar_axis_unit": xaxis_unit,
            "legend_description": " | ".join(extra_axes),
        },
    )
    apply_plot_settings(figure, plot_settings, theme)

    return figure


def _resolve_x_axis(x_axis: str | None, parameter_columns: list[str]) -> str:
    """Resolve a selector field name against namespaced result columns."""
    if x_axis in parameter_columns:
        return x_axis

    if x_axis:
        matching_columns = [column for column in parameter_columns if column.rsplit(":", 1)[-1] == x_axis]
        if matching_columns:
            return matching_columns[0]

    return parameter_columns[0]


def apply_plot_settings(
    figure: go.Figure,
    plot_settings: dict[str, Any] | None = None,
    theme: str = "light",
    polar_allowed: bool | None = None,
) -> go.Figure:
    """Apply persisted visual preferences to an existing Plotly figure."""
    settings = {**DEFAULT_PLOT_SETTINGS, **(plot_settings or {})}
    template = settings["template"]
    if template == "match-theme":
        template = "plotly_dark" if theme == "dark" else "plotly_white"
    is_dark = template == "plotly_dark"
    figure.update_layout(
        template=template,
        font={"size": settings["font_size"]},
        height=max(300, int(settings["graph_height"])),
        showlegend=bool(settings["show_legend"]),
        legend={
            "orientation": "h",
            "x": 0,
            "xanchor": "left",
            "y": -0.38,
            "yanchor": "top",
            "bgcolor": "rgba(32, 37, 43, 0.82)" if is_dark else "rgba(255, 255, 255, 0.92)",
            "bordercolor": "rgba(255, 255, 255, 0.24)" if is_dark else "rgba(96, 110, 123, 0.28)",
            "borderwidth": 1,
            "entrywidthmode": "pixels",
            "entrywidth": 170,
            "font": {"size": max(10, int(settings["font_size"] * 0.85)), "color": "#f3f5f7" if is_dark else "#26323b"},
        },
        margin={"b": 125},
        title_text="",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        plot_bgcolor="#20252b" if is_dark else "white",
        xaxis={"showgrid": bool(settings["show_grid"]), "type": settings["x_scale"]},
        yaxis={"showgrid": bool(settings["show_grid"]), "type": "log" if settings["log_y"] else "linear"},
    )
    metadata = figure.layout.meta if isinstance(figure.layout.meta, dict) else {}
    legend_description = metadata.get("legend_description")
    if legend_description and settings["show_legend"]:
        figure.add_annotation(
            text=legend_description,
            x=0,
            xref="paper",
            y=-0.29,
            yref="paper",
            xanchor="left",
            yanchor="top",
            showarrow=False,
            font={"size": max(10, int(settings["font_size"] * 0.78)), "color": "#f3f5f7" if is_dark else "#52616d"},
        )
    for trace in figure.data:
        if getattr(trace, "line", None) is not None:
            trace.line.width = settings["line_width"]

    if figure.layout.polar is not None:
        figure.update_layout(
            polar={
                "angularaxis": {"showgrid": bool(settings["show_grid"])},
                "radialaxis": {"showgrid": bool(settings["show_grid"]), "type": "log" if settings["log_y"] else "linear"},
            }
        )
    if polar_allowed is None:
        metadata = figure.layout.meta if isinstance(figure.layout.meta, dict) else {}
        polar_allowed = _is_angular_unit(metadata.get("polar_axis_unit"))

    if settings["coordinate_system"] == "polar" and polar_allowed:
        _convert_traces_to_polar(figure, bool(settings["show_grid"]), bool(settings["log_y"]))
    return figure


def _is_angular_unit(unit: Any) -> bool:
    """Return whether a serialized unit represents radians or degrees."""
    normalized = str(unit or "").strip().lower().replace("°", "degree")
    return normalized in {"rad", "radian", "radians", "deg", "degree", "degrees"}


def _polar_theta_degrees(values: Any, unit: Any) -> np.ndarray:
    """Convert serialized angular values to the degree convention used by Plotly."""
    theta = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    normalized = str(unit or "").strip().lower().replace("°", "degree")
    if normalized in {"rad", "radian", "radians"}:
        theta = np.rad2deg(theta)
    return np.mod(theta, 360.0)


def _format_plot_value(value: Any) -> str:
    """Keep legend values readable without exposing floating-point noise."""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)

    if not np.isfinite(numeric):
        return str(value)
    return f"{numeric:.4g}"


def _convert_traces_to_polar(figure: go.Figure, show_grid: bool, log_y: bool) -> None:
    """Convert compatible line traces to a polar subplot in place."""
    converted_traces = []
    converted = False
    for trace in figure.data:
        if trace.type not in {"scatter", "scattergl"} or trace.y is None:
            converted_traces.append(trace)
            continue

        # A Cartesian Scatter cannot be updated with ``r`` and ``theta``;
        # Plotly requires a Scatterpolar trace for those properties. Copy the
        # existing trace attributes so styling, hover text, and names survive.
        properties = trace.to_plotly_json()
        properties.pop("type", None)
        x_values = properties.pop("x", None)
        y_values = properties.pop("y")
        properties["r"] = list(y_values)
        properties["theta"] = list(x_values) if x_values is not None else list(
            np.linspace(0, 360, len(y_values), endpoint=False)
        )
        polar_properties = go.Scatterpolar._valid_props
        converted_traces.append(
            go.Scatterpolar(**{key: value for key, value in properties.items() if key in polar_properties})
        )
        converted = True

    if not converted:
        return

    # Plotly only permits assigning existing trace objects to ``figure.data``.
    # Clear the original traces, then add the newly constructed polar traces.
    figure.data = ()
    figure.add_traces(converted_traces)
    figure.update_layout(
        polar={
            "angularaxis": {"showgrid": show_grid},
            "radialaxis": {"showgrid": show_grid, "type": "log" if log_y else "linear"},
        },
        xaxis_visible=False,
        yaxis_visible=False,
    )


def build_summary(result: dict[str, Any] | None) -> list[dict[str, str]]:
    """Build compact summary metadata for the dashboard result cards."""
    if not result:
        return []

    LOGGER.debug("Building summary cards for measure=%s", result["measure"])

    return [
        {"label": "Measure", "value": result["measure"]},
        {"label": "Rows", "value": str(result["row_count"])},
        {"label": "Axes", "value": ", ".join(result["parameter_columns"]) or "none"},
    ]


def _parse_section_fields(field_specs: tuple[Any, ...], raw_values: Dict[str, Any]) -> Dict[str, Any]:
    """Parse raw text inputs according to the schema for one section."""
    parsed_values: Dict[str, Any] = {}

    for field in field_specs:
        raw_value = raw_values.get(field.name, field.default)

        LOGGER.debug("Parsing field %s raw_value=%r optional=%s", field.name, raw_value, field.optional)

        if field.optional and (raw_value is None or str(raw_value).strip() == ""):
            LOGGER.debug("Skipping optional empty field %s", field.name)
            continue

        parsed_value = _parse_field_value(field.kind, raw_value, field.unit)
        if field.name in _POSITIVE_FIELDS:
            _validate_positive_field(field.name, parsed_value)
        parsed_values[field.name] = parsed_value

    return parsed_values


def _validate_positive_field(name: str, value: Any) -> None:
    """Reject empty or non-positive values before they reach the C++ layer."""
    if value is None:
        raise ValueError(f"{name} must be positive.")

    magnitudes = getattr(value, "magnitude", value)
    values = np.asarray(magnitudes, dtype=float)
    if values.size == 0 or not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError(f"{name} must be positive.")


def _infer_variable_fields_for_section(field_specs: tuple[Any, ...], raw_values: Dict[str, Any]) -> list[str]:
    """Return field names that currently expand to multiple values in one section."""
    variable_fields: list[str] = []

    for field in field_specs:
        raw_value = raw_values.get(field.name, field.default)

        if field.optional and (raw_value is None or str(raw_value).strip() == ""):
            continue

        try:
            parsed_value = _parse_field_value(field.kind, raw_value, field.unit)
        except Exception:
            LOGGER.debug("Skipping variable-field inference for invalid field %s=%r", field.name, raw_value, exc_info=True)
            continue

        if _value_cardinality(parsed_value) > 1:
            variable_fields.append(field.name)

    return variable_fields


def _value_cardinality(value: Any) -> int:
    """Return how many effective values an already parsed field contains."""
    if value is None:
        return 0

    if isinstance(value, (str, bytes)):
        return 1

    try:
        return len(value)
    except TypeError:
        return 1


def _parse_field_value(kind: str, raw_value: Any, unit: Any) -> Any:
    """Dispatch parsing based on field type."""
    LOGGER.debug("Dispatching parser for kind=%s raw_value=%r unit=%r", kind, raw_value, unit)
    if kind == "quantity":
        return parse_quantity_expression(raw_value, unit)
    if kind == "numeric":
        return parse_numeric_expression(raw_value)
    if kind == "integer":
        return parse_numeric_expression(raw_value, integer=True)
    if kind == "polarization":
        return parse_polarization(raw_value, unit)
    if kind == "material":
        return parse_material_values(raw_value, medium=False)
    if kind == "medium":
        return parse_material_values(raw_value, medium=True)
    if kind == "mode":
        return parse_mode_numbers(raw_value)

    raise ValueError(f"Unsupported field kind '{kind}'.")
