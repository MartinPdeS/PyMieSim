"""Specialized dataframe utilities for PyMieSim experiment results."""

from typing import Self
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from MPSPlots.styles import scientific
from PyMieSim.material import ConstantMaterial, ConstantMedium


class PyMieSimDataFrame(pd.DataFrame):
    """Pandas dataframe subclass carrying PyMieSim units and plotting helpers.

    Instances behave like a regular :class:`pandas.DataFrame`, while storing
    unit metadata in ``attrs["units"]`` and exposing plotting utilities that
    understand PyMieSim material, medium and quantity objects.
    """

    @property
    def _constructor(self):
        """Return the dataframe subclass preserved by pandas operations.

        Returns
        -------
        type[PyMieSimDataFrame]
            The subclass constructor used by pandas internals.
        """
        return PyMieSimDataFrame

    # ---------------------------------------------------------
    # utilities
    # ---------------------------------------------------------

    def _units(self):
        """Return the mapping of dataframe columns to their registered units.

        Returns
        -------
        dict
            Mapping between column names and unit objects stored in
            ``attrs["units"]``.
        """
        return self.attrs.get("units", {})

    def _validate_column(self, column: str):
        """Raise a ``ValueError`` when ``column`` is not part of the dataframe.

        Parameters
        ----------
        column:
            Name of the dataframe column to validate.

        Returns
        -------
        None
            This method only validates input and raises on failure.
        """

        if column not in self.columns:
            available = ", ".join(self.columns)
            raise ValueError(
                f"Column '{column}' not found. Available columns: {available}"
            )

    def _measure_columns(self):
        """Return columns that look like measured outputs rather than parameters.

        Returns
        -------
        list[str]
            Unit-bearing columns that do not belong to scatterer, detector, or
            source parameter namespaces.
        """

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
        """Build an axis label using the stored unit metadata when available.

        Parameters
        ----------
        column:
            Column name used on one plot axis.

        Returns
        -------
        str
            Column label optionally suffixed with its stored unit.
        """

        units = self._units()

        if column in units:
            return f"{column} [{units[column]}]"

        return column

    def _clean_parameter_name(self, name):
        """Strip the namespace prefix from a column name.

        Parameters
        ----------
        name:
            Fully qualified PyMieSim column name.

        Returns
        -------
        str
            Final token after the ``:`` separator.
        """
        return name.split(":")[-1]

    def _format_parameter_label(self, parameters, values, precision=None):
        """Convert grouped parameter values into a compact legend label.

        Parameters
        ----------
        parameters:
            Iterable of grouped parameter column names.
        values:
            Group key values associated with ``parameters``.
        precision:
            Optional number of decimal places for numeric values.

        Returns
        -------
        str | None
            Joined legend label, or ``None`` when no parameters are provided.
        """

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

            if isinstance(parameter_value, (ConstantMaterial, ConstantMedium)):
                parameter_value = parameter_value.refractive_index

            elif hasattr(parameter_value, "name"):
                parameter_value = parameter_value.name

            elif hasattr(parameter_value, "magnitude"):
                parameter_value = parameter_value.magnitude

            if isinstance(parameter_value, complex):
                if np.isclose(parameter_value.imag, 0.0, atol=1e-15, rtol=0.0):
                    parameter_value = parameter_value.real

            if isinstance(parameter_value, (int, float, np.integer, np.floating)) and precision is not None:
                parameter_value = f"{parameter_value:.{precision}f}"

            parts.append(f"{clean_name}: {parameter_value}")

        return " | ".join(parts)

    # ---------------------------------------------------------
    # unit scaling utilities
    # ---------------------------------------------------------

    def scale_column_unit(self, column: str, target_unit) -> Self:
        """
        Return a copy with one column converted to ``target_unit``.

        Parameters
        ----------
        column:
            Name of the column to convert.
        target_unit:
            Pint-compatible unit used as the destination unit.

        Returns
        -------
        PyMieSimDataFrame
            A copy of the dataframe with converted values and updated unit
            metadata.

        Raises
        ------
        ValueError
            If the column is missing or has no registered unit.
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
        Return a copy with each unit-bearing column converted to a compact unit.

        Parameters
        ----------
        None

        The compact unit is chosen from the largest absolute finite value found
        in each column, which keeps plots and tables easier to read without
        changing the represented physical quantity.

        Returns
        -------
        PyMieSimDataFrame
            Copy of the dataframe with compacted display units.
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

    def _resolve_projection(self, projection=None, project=None):
        """Normalize accepted projection names and the legacy ``project`` alias.

        Parameters
        ----------
        projection:
            Explicit projection name passed by the caller.
        project:
            Legacy alias for ``projection`` kept for backward compatibility.

        Returns
        -------
        str | None
            ``"polar"`` for polar axes, otherwise ``None`` for regular axes.
        """

        resolved_projection = projection if projection is not None else project

        if resolved_projection in (None, "normal", "rectilinear", "cartesian"):
            return None

        if resolved_projection != "polar":
            raise ValueError(
                "projection must be one of None, 'normal', 'rectilinear', "
                "'cartesian', or 'polar'."
            )

        return "polar"

    def _create_figure(self, axes=None, projection=None, figure_size=None):
        """Create a new figure or reuse an existing axes object for plotting.

        Parameters
        ----------
        axes:
            Existing matplotlib axes to reuse.
        projection:
            Projection name returned by :meth:`_resolve_projection`.
        figure_size:
            Requested figure size in inches.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Figure and axes used for plotting.
        """

        if axes is not None:
            if projection == "polar" and getattr(axes, "name", None) != "polar":
                raise ValueError(
                    "A polar projection was requested but the provided axes is not polar."
                )

            return axes.figure, axes

        subplot_kwargs = {"projection": projection} if projection else None

        return plt.subplots(1, 1, figsize=figure_size, subplot_kw=subplot_kwargs)

    def _finalize_figure(
        self,
        figure,
        title=None,
        tight_layout=True,
        save_as=None,
        show=True,
        xscale=None,
        yscale=None,
        xlim=None,
        ylim=None,
    ):
        """Apply final figure-wide styling, limits, layout, saving, and display.

        Parameters
        ----------
        figure:
            Matplotlib figure produced by one of the plotting methods.
        title:
            Optional title applied to the single axes or the whole figure.
        tight_layout:
            Whether to call ``figure.tight_layout()``.
        save_as:
            Optional output path used to save the figure.
        show:
            Whether to display the figure with ``plt.show()``.
        xscale, yscale:
            Optional Matplotlib axis scales.
        xlim, ylim:
            Optional axis limits.

        Returns
        -------
        matplotlib.figure.Figure
            Finalized figure ready for display or further customization.
        """

        for ax in figure.axes:
            if xscale is not None:
                ax.set_xscale(xscale)

            if yscale is not None:
                ax.set_yscale(yscale)

            if xlim is not None:
                ax.set_xlim(*xlim)

            if ylim is not None:
                ax.set_ylim(*ylim)

        if title:
            if len(figure.axes) == 1:
                figure.axes[0].set_title(title)
            else:
                figure.suptitle(title)

        if tight_layout:
            figure.tight_layout()

        if save_as is not None:
            figure.savefig(save_as, dpi=300)

        if show:
            plt.show()

        return figure

    def _apply_axis_resolution(self, ax, x_precision, y_precision):
        """Format numeric axis ticks with a fixed decimal precision.

        Parameters
        ----------
        ax:
            Matplotlib axes whose tick formatter should be updated.
        x_precision:
            Number of decimals for the x-axis ticks, or ``None`` to keep the
            default formatter.
        y_precision:
            Number of decimals for the y-axis ticks, or ``None`` to keep the
            default formatter.

        Returns
        -------
        None
            Tick formatters are updated in place.
        """

        if x_precision is not None:
            ax.xaxis.set_major_formatter(
                FormatStrFormatter(f"%.{x_precision}f")
            )

        if y_precision is not None:
            ax.yaxis.set_major_formatter(
                FormatStrFormatter(f"%.{y_precision}f")
            )

    def _format_legend(self, ax, precision):
        """Reformat numeric values in legend labels with a shared precision.

        Parameters
        ----------
        ax:
            Matplotlib axes containing the legend to update.
        precision:
            Number of decimals used for numeric values inside the legend.

        Returns
        -------
        None
            Legend text objects are updated in place when present.
        """

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
        """Plot one or more measured columns against an experiment parameter.

        Parameters
        ----------
        *args:
            Positional arguments forwarded to :meth:`plot_standard` or
            :meth:`plot_with_std`.
        **kwargs:
            Keyword arguments forwarded to :meth:`plot_standard` or
            :meth:`plot_with_std`.

        When ``std`` is provided the dataframe is first aggregated across
        repeated samples and the mean curve is plotted with a shaded standard
        deviation envelope. Otherwise each curve is plotted directly from the
        dataframe values.

        This dispatcher forwards all arguments to :meth:`plot_standard` or
        :meth:`plot_with_std`.

        Returns
        -------
        matplotlib.figure.Figure
            Figure returned by the selected plotting backend.
        """

        if kwargs.get("std", None) is not None:
            return self.plot_with_std(*args, **kwargs)

        return self.plot_standard(*args, **kwargs)

    # ---------------------------------------------------------
    # standard plotting
    # ---------------------------------------------------------
    def plot_standard(
        self,
        x: str,
        y: list[str] | None = None,
        legend_resolution: int | None = 1,
        x_tick_resolution: int | None = 1,
        y_tick_resolution: int | None = 1,
        title: str | None = None,
        figure_size: tuple[float, float] = (10, 6),
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        xscale: str | None = None,
        yscale: str | None = None,
        projection: str | None = None,
        project: str | None = None,
        axes=None,
        tight_layout: bool = True,
        save_as: str | None = None,
        show: bool = True,
        style=scientific,
        **kwargs,
    ):
        """Plot raw dataframe curves without aggregation.

        Parameters
        ----------
        x:
            Column used for the x-axis.
        y:
            Measured columns to plot. When omitted, every measure-like column is
            plotted.
        legend_resolution:
            Number of decimals used when formatting numeric values in legend
            entries.
        x_tick_resolution, y_tick_resolution:
            Number of decimals displayed on numeric axis tick labels.
        title:
            Optional plot title.
        figure_size:
            Figure size in inches as ``(width, height)``.
        xlim, ylim:
            Optional axis limits.
        xscale, yscale:
            Matplotlib scale names such as ``"linear"`` or ``"log"``.
        projection:
            Plot projection. Use ``"polar"`` for polar axes or keep ``None`` for
            a standard Cartesian plot.
        project:
            Backward-compatible alias for ``projection``.
        axes:
            Existing Matplotlib axes to draw on. When omitted, a new figure and
            axes pair is created.
        tight_layout:
            Whether to call ``figure.tight_layout()`` before returning.
        save_as:
            Optional file path used to save the figure.
        show:
            Whether to call :func:`matplotlib.pyplot.show`.
        style:
            Matplotlib style applied during figure creation. Defaults to
            ``MPSPlots.styles.scientific``.
        **kwargs:
            Additional keyword arguments forwarded to ``Axes.plot``.

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import warnings

        resolved_projection = self._resolve_projection(
            projection=projection,
            project=project,
        )

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
            """Convert PyMieSim objects and quantities into scalar plot values.

            Parameters
            ----------
            value:
                Raw dataframe entry to normalize for plotting.

            Returns
            -------
            object
                Scalar-compatible value extracted from the original object.
            """
            if isinstance(value, (ConstantMaterial, ConstantMedium)):
                return value.refractive_index

            if hasattr(value, "name"):
                return value.name

            if hasattr(value, "magnitude"):
                return value.magnitude

            return value

        def _normalize_numeric_series(series, use_abs_for_complex=False, warn_on_complex_real_projection=False):
            """Convert a series into floats while handling complex-valued entries.

            Parameters
            ----------
            series:
                Pandas series to convert into numeric plot coordinates.
            use_abs_for_complex:
                Whether to convert complex values to their magnitude.
            warn_on_complex_real_projection:
                Whether to keep track of non-zero imaginary components when the
                real part is used for plotting.

            Returns
            -------
            tuple[numpy.ndarray | None, bool]
                Numeric array when conversion succeeds and a flag indicating
                whether a discarded imaginary component was detected.
            """
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

        def _normalize_categorical_series(series):
            """Convert a series to string labels suitable for categorical plotting.

            Parameters
            ----------
            series:
                Pandas series containing categorical-like values.

            Returns
            -------
            pandas.Series
                String-valued series suitable for axis tick labels.
            """
            return series.map(lambda value: str(_extract_scalar(value)))

        def _compose_curve_label(parameter_label, measure_name):
            """Build the legend label for a plotted measure curve.

            Parameters
            ----------
            parameter_label:
                Legend fragment describing the grouped parameter values.
            measure_name:
                Name of the measure currently being plotted.

            Returns
            -------
            str
                Full legend label for a single plotted line.
            """
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
            categorical_labels = _normalize_categorical_series(df[x])
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

        with plt.style.context(style):
            figure, ax = self._create_figure(
                axes=axes,
                projection=resolved_projection,
                figure_size=figure_size,
            )

            groups = list(df.groupby(parameters, sort=False, dropna=False)) if parameters else [(None, df)]

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
                    group_labels = _normalize_categorical_series(group[x])
                    x_values = group_labels.map(category_to_position).to_numpy(dtype=float)

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

            if x_is_numeric:
                self._apply_axis_resolution(ax, x_tick_resolution, y_tick_resolution)
            else:
                self._apply_axis_resolution(ax, None, y_tick_resolution)

            self._format_legend(ax, legend_resolution)

            return self._finalize_figure(
                figure=figure,
                title=title,
                tight_layout=tight_layout,
                save_as=save_as,
                show=show,
                xscale=xscale,
                yscale=yscale,
                xlim=xlim,
                ylim=ylim,
            )
    # ---------------------------------------------------------
    # std plotting
    # ---------------------------------------------------------

    def plot_with_std(
        self,
        x: str,
        std: str,
        y: list[str] | None = None,
        alpha: float = 0.4,
        legend_resolution: int = 4,
        x_tick_resolution: int = 4,
        y_tick_resolution: int = 4,
        title: str | None = None,
        figure_size: tuple[float, float] | None = None,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        xscale: str | None = None,
        yscale: str | None = None,
        projection: str | None = None,
        project: str | None = None,
        axes=None,
        tight_layout: bool = True,
        save_as: str | None = None,
        show: bool = True,
        style=scientific,
        **kwargs,
    ):
        """Plot mean curves with a shaded standard deviation envelope.

        Parameters
        ----------
        x:
            Column used for the x-axis.
        std:
            Column that defines replicates before aggregation. Data are grouped
            by ``x`` and every varying non-measure column except ``std``.
        y:
            Measured columns to aggregate and plot. When omitted, every
            measure-like column is included.
        alpha:
            Opacity of the shaded standard deviation band.
        legend_resolution:
            Number of decimals used when formatting numeric values in legend
            entries.
        x_tick_resolution, y_tick_resolution:
            Number of decimals displayed on numeric axis tick labels.
        title, figure_size, xlim, ylim, xscale, yscale, projection, project,
        axes, tight_layout, save_as, show, style:
            Same behavior as in :meth:`plot_standard`.
        **kwargs:
            Additional keyword arguments forwarded to ``Axes.plot``.

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import warnings

        resolved_projection = self._resolve_projection(
            projection=projection,
            project=project,
        )

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
            """Convert PyMieSim objects and quantities into scalar plot values.

            Parameters
            ----------
            value:
                Raw dataframe entry to normalize for plotting.

            Returns
            -------
            object
                Scalar-compatible value extracted from the original object.
            """
            if isinstance(value, (ConstantMaterial, ConstantMedium)):
                return value.refractive_index

            if hasattr(value, "name"):
                return value.name

            if hasattr(value, "magnitude"):
                return value.magnitude

            return value

        def _to_real_array(series, use_abs_for_complex=False, warn_on_complex_real_projection=False):
            """Convert a series into floats while handling complex-valued entries.

            Parameters
            ----------
            series:
                Pandas series to convert into numeric plot coordinates.
            use_abs_for_complex:
                Whether to convert complex values to their magnitude.
            warn_on_complex_real_projection:
                Whether to keep track of non-zero imaginary components when the
                real part is used for plotting.

            Returns
            -------
            tuple[numpy.ndarray | None, bool]
                Numeric array when conversion succeeds and a flag indicating
                whether a discarded imaginary component was detected.
            """
            numeric_values = []
            has_nonzero_imaginary_part = False

            for value in series.to_numpy():
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
                else:
                    try:
                        numeric_values.append(float(scalar_value))
                    except (TypeError, ValueError):
                        return None, has_nonzero_imaginary_part

            return np.asarray(numeric_values, dtype=float), has_nonzero_imaginary_part

        def _to_categorical_series(series):
            """Convert a series to string labels suitable for categorical plotting.

            Parameters
            ----------
            series:
                Pandas series containing categorical-like values.

            Returns
            -------
            pandas.Series
                String-valued series suitable for axis tick labels.
            """
            return series.map(lambda value: str(_extract_scalar(value)))

        def _compose_curve_label(parameter_label, measure_name):
            """Build the legend label for a plotted mean curve.

            Parameters
            ----------
            parameter_label:
                Legend fragment describing the grouped parameter values.
            measure_name:
                Name of the aggregated measure currently being plotted.

            Returns
            -------
            str
                Full legend label for a single plotted line.
            """
            clean_measure_name = self._clean_parameter_name(measure_name)

            if len(y) == 1:
                return parameter_label if parameter_label else clean_measure_name

            if parameter_label:
                return f"{clean_measure_name} | {parameter_label}"

            return clean_measure_name

        x_numeric_values, x_has_nonzero_imaginary_part = _to_real_array(
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
            categorical_labels = _to_categorical_series(df[x])
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

            std_values, _ = _to_real_array(group_df[std], use_abs_for_complex=False)
            row[f"{std}__mean"] = std_values.mean()
            row[f"{std}__std"] = std_values.std(ddof=1) if len(group_df) > 1 else 0.0

            for measure in y:
                values, _ = _to_real_array(group_df[measure], use_abs_for_complex=True)
                row[f"{measure}__mean"] = values.mean()
                row[f"{measure}__std"] = values.std(ddof=1) if len(values) > 1 else 0.0

            aggregated_rows.append(row)

        aggregated_df = pd.DataFrame(aggregated_rows)

        with plt.style.context(style):
            figure, ax = self._create_figure(
                axes=axes,
                projection=resolved_projection,
                figure_size=figure_size,
            )

            grouped_curves = list(aggregated_df.groupby(parameters, sort=False, dropna=False)) if parameters else [(None, aggregated_df)]

            for key, group_df in grouped_curves:
                if x_is_numeric:
                    x_values, _ = _to_real_array(
                        group_df[x],
                        warn_on_complex_real_projection=True,
                    )
                    sorting_index = np.argsort(x_values)

                    group_df = group_df.iloc[sorting_index]
                    x_values = x_values[sorting_index]
                else:
                    group_labels = _to_categorical_series(group_df[x])
                    x_values = group_labels.map(category_to_position).to_numpy(dtype=float)

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

            if x_is_numeric:
                self._apply_axis_resolution(ax, x_tick_resolution, y_tick_resolution)
            else:
                self._apply_axis_resolution(ax, None, y_tick_resolution)

            self._format_legend(ax, legend_resolution)

            return self._finalize_figure(
                figure=figure,
                title=title,
                tight_layout=tight_layout,
                save_as=save_as,
                show=show,
                xscale=xscale,
                yscale=yscale,
                xlim=xlim,
                ylim=ylim,
            )

    def get(self, key: str):
        """Return a column as a NumPy array multiplied by its registered unit.

        Parameters
        ----------
        key:
            Name of the dataframe column to extract.

        Returns
        -------
        numpy.ndarray | pint.Quantity
            Column values multiplied by their registered unit when available.
        """
        return self[key].to_numpy() * self._units().get(key, 1)

    # ---------------------------------------------------------
    # dataframe utilities
    # ---------------------------------------------------------

    def add_weight(self, weight_column: str, weight: np.ndarray) -> Self:
        """Return a copy with ``weight_column`` scaled by ``weight`` elementwise.

        Parameters
        ----------
        weight_column:
            Name of the column to rescale.
        weight:
            Array of multiplicative weights broadcastable to the column shape.

        Returns
        -------
        PyMieSimDataFrame
            Copy of the dataframe with weighted values.
        """

        self._validate_column(weight_column)

        df = self.copy()

        df[weight_column] = df[weight_column] * weight.squeeze()

        return df

    def sum_over(self, column: str) -> Self:
        """Group rows by ``column`` and return the summed dataframe.

        Parameters
        ----------
        column:
            Name of the grouping column.

        Returns
        -------
        PyMieSimDataFrame
            Summed dataframe indexed by ``column``.
        """

        self._validate_column(column)

        return self.groupby(column).sum()