from typing import Self
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from MPSPlots import helper
from PyMieSim.units import ureg


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

        parts = []

        for p, v in zip(parameters, values):

            p = self._clean_parameter_name(p)

            if isinstance(v, (int, float)) and precision is not None:
                v = f"{v:.{precision}f}"

            parts.append(f"{p}: {v}")

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

            if values.size == 0:
                continue

            max_val = np.nanmax(np.abs(values))

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

                    name, value = item.split(":")
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
        y = [y] if isinstance(y, str) else y
        self._validate_column(x)

        df = self.copy()

        measures = self._measure_columns()

        if y is None:
            y = measures

        parameters = [
            c for c in df.columns
            if c not in ([x] + y)
            and not df[c].isna().all()
            and df[c].nunique(dropna=True) > 1
        ]

        figure, ax = self._create_figure()

        groups = df.groupby(parameters) if parameters else [(None, df)]

        for key, group in groups:
            group = group.sort_values(x)

            label = self._format_parameter_label(
                parameters,
                key,
                precision=legend_resolution
            )

            for measure in y:
                ax.plot(
                    group[x].to_numpy(),
                    group[measure].to_numpy(),
                    label=label if label else measure,
                    **kwargs,
                )

        ax.set_xlabel(self._axis_label(x))
        ax.set_ylabel(self._axis_label(y[0]))

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

        parameters = [
            c for c in df.columns
            if c not in ([x, std] + y)
        ]

        grouped = df.groupby(parameters + [x])

        mean = grouped.mean().reset_index()
        stdv = grouped.std().reset_index()

        figure, ax = self._create_figure()

        groups = mean.groupby(parameters) if parameters else [(None, mean)]

        for measure in y:

            for key, group_df in groups:

                group_df = group_df.sort_values(x)

                if parameters:
                    mask = stdv[parameters].eq(group_df[parameters].iloc[0]).all(axis=1)
                    std_group = stdv[mask]
                else:
                    std_group = stdv

                label = self._format_parameter_label(
                    parameters,
                    key,
                    precision=legend_resolution
                )

                ax.plot(
                    group_df[x],
                    group_df[measure],
                    label=label if label else measure,
                    **kwargs,
                )

                ax.fill_between(
                    group_df[x],
                    group_df[measure] - std_group[measure],
                    group_df[measure] + std_group[measure],
                    alpha=alpha,
                )

        ax.set_xlabel(self._axis_label(x))
        ax.set_ylabel(self._axis_label(y[0]))

        ax.legend()

        self._apply_axis_resolution(ax, x_tick_resolution, y_tick_resolution)
        self._format_legend(ax, legend_resolution)

        return figure

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