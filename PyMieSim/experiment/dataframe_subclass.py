from typing import Self
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from MPSPlots import helper
from PyMieSim.material import ConstantMaterial, ConstantMedium


class PyMieSimDataFrame(pd.DataFrame):

    @property
    def _constructor(self):
        return PyMieSimDataFrame

    # ---------------------------------------------------------
    # utilities
    # ---------------------------------------------------------

    def _units(self):
        return self.attrs.get("units", {})

    def _validate_column(self, column: str):

        if column not in self.columns:
            available = ", ".join(self.columns)
            raise ValueError(
                f"Column '{column}' not found. Available columns: {available}"
            )

    def _measure_columns(self):

        units = self._units()

        measures = []

        for column in self.columns:

            if column in units:

                if not (
                    column.startswith("scatterer:")
                    or column.startswith("detector:")
                    or column.startswith("source:")
                ):
                    measures.append(column)

        return measures

    def _axis_label(self, column):

        units = self._units()

        if column in units:
            return f"{column} [{units[column]}]"

        return column

    def _clean_parameter_name(self, name):
        return name.split(":")[-1]

    def _format_parameter_label(self, parameters, values, precision=None):

        if parameters is None:
            return None

        if not isinstance(values, tuple):
            values = (values,)

        units = self._units()
        parts = []

        for parameter_name, parameter_value in zip(parameters, values):

            clean_name = self._clean_parameter_name(parameter_name)

            unit = units.get(parameter_name, None)
            if unit is not None:
                clean_name = f"{clean_name} [{unit:~}]"

            if isinstance(parameter_value, (int, float, np.integer, np.floating)) and precision is not None:
                parameter_value = f"{parameter_value:.{precision}f}"

            parts.append(f"{clean_name}: {parameter_value}")

        return " | ".join(parts)

    # ---------------------------------------------------------
    # unit scaling utilities
    # ---------------------------------------------------------

    def scale_column_unit(self, column: str, target_unit) -> Self:
        """
        Convert a column to a new unit.
        """

        self._validate_column(column)

        units = self._units()

        if column not in units:
            raise ValueError(f"No unit registered for column '{column}'")

        current_unit = units[column]

        factor = (1 * current_unit).to(target_unit).magnitude

        df = self.copy()

        df[column] = df[column] * factor

        df.attrs["units"][column] = target_unit

        return df

    def to_compact(self) -> Self:
        """
        Scale each column to its most compact unit based on the largest value.
        """

        df = self.copy()

        units = df._units()

        for column, unit in units.items():

            values = df[column].to_numpy()

            try:
                max_val = np.nanmax(np.abs(values))
            except:
                continue

            if max_val == 0:
                continue

            quantity = max_val * unit

            compact_unit = quantity.to_compact().units

            factor = (1 * unit).to(compact_unit).magnitude

            df[column] = df[column] * factor

            df.attrs["units"][column] = compact_unit

        return df

    # ---------------------------------------------------------
    # plotting helpers
    # ---------------------------------------------------------

    def _create_figure(self):
        return plt.subplots(1, 1)

    def _apply_axis_resolution(self, ax, x_precision, y_precision):

        if x_precision is not None:
            ax.xaxis.set_major_formatter(
                FormatStrFormatter(f"%.{x_precision}f")
            )

        if y_precision is not None:
            ax.yaxis.set_major_formatter(
                FormatStrFormatter(f"%.{y_precision}f")
            )

    def _format_legend(self, ax, precision):

        legend = ax.get_legend()

        if legend is None or precision is None:
            return

        for text in legend.get_texts():

            label = text.get_text()
            parts = []

            for item in label.split(" | "):

                if ":" in item:

                    name, value = item.split(":", 1)
                    name = name.strip()
                    value = value.strip()

                    try:
                        value = float(value)
                        value = f"{value:.{precision}f}"
                        parts.append(f"{name}: {value}")
                    except ValueError:
                        parts.append(item)

                else:
                    parts.append(item)

            text.set_text(" | ".join(parts))

    # ---------------------------------------------------------
    # plotting dispatcher
    # ---------------------------------------------------------

    def plot(self, *args, **kwargs):

        if kwargs.get("std", None) is not None:
            return self.plot_with_std(*args, **kwargs)

        return self.plot_standard(*args, **kwargs)

    # ---------------------------------------------------------
    # standard plotting
    # ---------------------------------------------------------
    @helper.post_mpl_plot
    def plot_standard(
        self,
        x: str,
        y: list[str] | None = None,
        legend_resolution: int | None = 1,
        x_tick_resolution: int | None = 1,
        y_tick_resolution: int | None = 1,
        **kwargs,
    ):
        import warnings

        y = [y] if isinstance(y, str) else y
        self._validate_column(x)

        df = self.copy()

        measures = self._measure_columns()

        if y is None:
            y = measures

        for measure in y:
            self._validate_column(measure)

        parameters = [
            column for column in df.columns
            if column not in ([x] + y)
            and not df[column].isna().all()
            and df[column].nunique(dropna=True) > 1
        ]

        def _extract_scalar(value):
            if hasattr(value, "magnitude"):
                return value.magnitude
            return value

        def _normalize_numeric_series(series, use_abs_for_complex=False, warn_on_complex_real_projection=False):
            raw_values = series.to_numpy()
            numeric_values = []
            has_nonzero_imaginary_part = False

            for value in raw_values:
                scalar_value = _extract_scalar(value)

                if isinstance(scalar_value, complex) or np.iscomplexobj(scalar_value):
                    if warn_on_complex_real_projection:
                        if not np.isclose(np.imag(scalar_value), 0.0, atol=1e-15, rtol=0.0):
                            has_nonzero_imaginary_part = True
                        numeric_values.append(float(np.real(scalar_value)))
                    elif use_abs_for_complex:
                        numeric_values.append(float(np.abs(scalar_value)))
                    else:
                        numeric_values.append(float(np.real(scalar_value)))
                    continue

                try:
                    numeric_values.append(float(scalar_value))
                except (TypeError, ValueError):
                    return None, has_nonzero_imaginary_part

            return np.asarray(numeric_values, dtype=float), has_nonzero_imaginary_part

        def _compose_curve_label(parameter_label, measure_name):
            clean_measure_name = self._clean_parameter_name(measure_name)

            if len(y) == 1:
                return parameter_label if parameter_label else clean_measure_name

            if parameter_label:
                return f"{clean_measure_name} | {parameter_label}"

            return clean_measure_name

        x_numeric_values, x_has_nonzero_imaginary_part = _normalize_numeric_series(
            df[x],
            warn_on_complex_real_projection=True,
        )

        x_is_numeric = x_numeric_values is not None
        x_used_real_part = False

        if x_is_numeric:
            raw_x_values = df[x].to_numpy()
            x_used_real_part = any(
                isinstance(_extract_scalar(value), complex) or np.iscomplexobj(_extract_scalar(value))
                for value in raw_x_values
            )
        else:
            categorical_labels = df[x].astype(str)
            unique_labels = pd.Index(categorical_labels).unique()
            category_to_position = {
                label: float(index)
                for index, label in enumerate(unique_labels)
            }

        if x_has_nonzero_imaginary_part:
            warnings.warn(
                f"Column {x!r} contains complex values on the x axis with nonzero imaginary part. "
                "Using the real part for plotting.",
                UserWarning,
                stacklevel=3,
            )

        figure, ax = self._create_figure()

        groups = list(df.groupby(parameters, sort=False)) if parameters else [(None, df)]

        for key, group in groups:

            if x_is_numeric:
                group_x_values, _ = _normalize_numeric_series(
                    group[x],
                    warn_on_complex_real_projection=True,
                )
                sorting_index = np.argsort(group_x_values)
                group = group.iloc[sorting_index]
                x_values = group_x_values[sorting_index]
            else:
                x_values = group[x].astype(str).map(category_to_position).to_numpy(dtype=float)

            parameter_label = self._format_parameter_label(
                parameters,
                key,
                precision=legend_resolution,
            )

            for measure in y:
                y_values, _ = _normalize_numeric_series(
                    group[measure],
                    use_abs_for_complex=True,
                )

                if y_values is None:
                    raise TypeError(
                        f"Column {measure!r} could not be converted to real values for plotting."
                    )

                curve_label = _compose_curve_label(parameter_label, measure)

                ax.plot(
                    x_values,
                    y_values,
                    label=curve_label,
                    **kwargs,
                )

        if not x_is_numeric:
            ax.set_xticks(list(category_to_position.values()))
            ax.set_xticklabels(list(category_to_position.keys()), rotation=45, ha='right')

        x_label = self._axis_label(x)

        if len(y) == 1:
            y_label = self._axis_label(y[0])
        else:
            y_label = "Measure"

        if x_is_numeric and x_used_real_part:
            x_label = f"Re({x_label})"

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        ax.legend()

        self._apply_axis_resolution(ax, x_tick_resolution, y_tick_resolution)
        self._format_legend(ax, legend_resolution)

        return figure
    # ---------------------------------------------------------
    # std plotting
    # ---------------------------------------------------------

    @helper.post_mpl_plot
    def plot_with_std(
        self,
        x: str,
        std: str,
        y: list[str] | None = None,
        alpha: float = 0.4,
        legend_resolution: int = 4,
        x_tick_resolution: int = 4,
        y_tick_resolution: int = 4,
        **kwargs,
    ):
        self._validate_column(x)
        self._validate_column(std)

        df = self.copy()

        measures = self._measure_columns()

        if y is None:
            y = measures
        elif isinstance(y, str):
            y = [y]

        for measure in y:
            self._validate_column(measure)

        parameters = [
            column for column in df.columns
            if column not in ([x, std] + y)
            and not df[column].isna().all()
            and df[column].nunique(dropna=True) > 1
        ]

        def _extract_scalar(value):
            if hasattr(value, "magnitude"):
                return value.magnitude
            return value

        def _to_real_array(series, use_abs_for_complex=False):
            numeric_values = []

            for value in series.to_numpy():
                scalar_value = _extract_scalar(value)

                if isinstance(scalar_value, complex) or np.iscomplexobj(scalar_value):
                    if use_abs_for_complex:
                        numeric_values.append(float(np.abs(scalar_value)))
                    else:
                        numeric_values.append(float(np.real(scalar_value)))
                else:
                    numeric_values.append(float(scalar_value))

            return np.asarray(numeric_values, dtype=float)

        def _compose_curve_label(parameter_label, measure_name):
            clean_measure_name = self._clean_parameter_name(measure_name)

            if len(y) == 1:
                return parameter_label if parameter_label else clean_measure_name

            if parameter_label:
                return f"{clean_measure_name} | {parameter_label}"

            return clean_measure_name

        aggregated_rows = []

        grouped = df.groupby(parameters + [x], sort=False, dropna=False)

        for group_key, group_df in grouped:
            if not isinstance(group_key, tuple):
                group_key = (group_key,)

            row = {}
            key_names = parameters + [x]

            for key_name, key_value in zip(key_names, group_key):
                row[key_name] = key_value

            row["_replicate_count"] = len(group_df)
            row[f"{std}__mean"] = _to_real_array(group_df[std]).mean()
            row[f"{std}__std"] = _to_real_array(group_df[std]).std(ddof=1) if len(group_df) > 1 else 0.0

            for measure in y:
                values = _to_real_array(group_df[measure], use_abs_for_complex=True)
                row[f"{measure}__mean"] = values.mean()
                row[f"{measure}__std"] = values.std(ddof=1) if len(values) > 1 else 0.0

            aggregated_rows.append(row)

        aggregated_df = pd.DataFrame(aggregated_rows)

        figure, ax = self._create_figure()

        grouped_curves = list(aggregated_df.groupby(parameters, sort=False, dropna=False)) if parameters else [(None, aggregated_df)]

        for key, group_df in grouped_curves:
            x_values = _to_real_array(group_df[x])
            sorting_index = np.argsort(x_values)

            group_df = group_df.iloc[sorting_index]
            x_values = x_values[sorting_index]

            parameter_label = self._format_parameter_label(
                parameters,
                key,
                precision=legend_resolution,
            )

            for measure in y:
                y_mean = group_df[f"{measure}__mean"].to_numpy(dtype=float)
                y_std = group_df[f"{measure}__std"].to_numpy(dtype=float)

                curve_label = _compose_curve_label(parameter_label, measure)

                ax.plot(
                    x_values,
                    y_mean,
                    label=curve_label,
                    **kwargs,
                )

                ax.fill_between(
                    x_values,
                    y_mean - y_std,
                    y_mean + y_std,
                    alpha=alpha,
                )

        ax.set_xlabel(self._axis_label(x))

        if len(y) == 1:
            ax.set_ylabel(self._axis_label(y[0]))
        else:
            ax.set_ylabel("Measure")

        ax.legend()

        self._apply_axis_resolution(ax, x_tick_resolution, y_tick_resolution)
        self._format_legend(ax, legend_resolution)

        return figure

    def get(self, key: str):
        return self[key].to_numpy() * self._units().get(key, 1)

    # ---------------------------------------------------------
    # dataframe utilities
    # ---------------------------------------------------------

    def add_weight(self, weight_column: str, weight: np.ndarray) -> Self:

        self._validate_column(weight_column)

        df = self.copy()

        df[weight_column] = df[weight_column] * weight.squeeze()

        return df

    def sum_over(self, column: str) -> Self:

        self._validate_column(column)

        return self.groupby(column).sum()