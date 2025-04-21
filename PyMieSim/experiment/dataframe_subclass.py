from MPSPlots.styles import mps
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional

class PyMieSimDataFrame(pd.DataFrame):
    """
    A subclass of pandas DataFrame with custom plot methods tailored for
    multi-indexed data containing pint-quantified values.
    """

    @property
    def _constructor(self):
        """
        Ensures that operations returning DataFrames return instances of PyMieSimDataFrame.
        """
        return PyMieSimDataFrame

    def _validate_axis(self, axis: str) -> None:
        # Validate that 'x' is a valid index level.
        if axis not in self.index.names:
            available = ", ".join(self.index.names)
            raise ValueError(f"x parameter '{axis}' is not in the DataFrame index. Available index levels: {available}")

    def _get_complementary_axis(self, *axis) -> tuple:
        return [level for level in self.index.names if level not in axis]

    def plot(self,
        x: str,
        std: Optional[str] = None,
        alpha: float = 0.4,
        ax: Optional[plt.Axes] = None,
        show: bool = True,
        xscale: bool = 'linear',
        yscale: bool = 'linear',
        **kwargs) -> plt.Axes:
        """
        Plots the DataFrame using a specified MultiIndex level for the x-axis.
        Optionally, if a standard deviation level is provided, it will also
        plot the mean ± std.

        Parameters
        ----------
        x : str
            The MultiIndex level to use for the x-axis.
        alpha : float, optional
            The transparency level for the shaded standard deviation region.
        std : Optional[str], optional
            The MultiIndex level used for standard deviation calculation.
        ax : Optional[plt.Axes], optional
            The matplotlib Axes on which to plot. If None, a new Axes is created.
        show : bool, optional
            If True, displays the plot.
        **kwargs : dict
            Additional keyword arguments passed to the underlying plotting functions.

        Returns
        -------
        plt.Axes
            The matplotlib Axes with the plot.
        """
        self._validate_axis(axis=x)

        if std is not None:
            self._validate_axis(axis=std)


        with plt.style.context(mps):
            if ax is None:
                _, ax = plt.subplots()

            ax.set(xscale=xscale, yscale=yscale)

            if std is not None:
                self._plot_with_std(ax=ax, x=x, std=std, alpha=alpha, **kwargs)
            else:
                self._plot_without_std(ax=ax, x=x, **kwargs)

            self._format_legend(ax)

            if show:
                plt.show()

            return ax

    def _plot_with_std(self, ax: plt.Axes, x: str, std: str, alpha: float = 0.5, **kwargs) -> None:
        """
        Plot the mean with standard deviation shading.
        Expects the DataFrame to have a MultiIndex and a 'pint' attribute.

        Parameters
        ----------
        ax : plt.Axes
            The matplotlib Axes to plot on.
        x : str
            The MultiIndex level to use for the x-axis.
        std : str
            The MultiIndex level used for standard deviation calculation.
        alpha : float, optional
            Transparency for the std shading.
        show : bool, optional
            Whether to call plt.show() after plotting.
        **kwargs : dict
            Additional keyword arguments for line styling.
        """

        # Determine levels to unstack
        no_std_levels = self._get_complementary_axis(std)
        no_x_levels = self._get_complementary_axis(x, std)

        # Calculate standard deviation and mean, preserving pint units
        std_df = (
            self.pint.dequantify()
            .unstack(no_std_levels)
            .apply(np.std)
            .to_frame(name="std")
            .unstack("unit")
            .pint.quantify()
        )
        mean_df = (
            self.pint.dequantify()
            .unstack(no_std_levels)
            .apply(np.mean)
            .to_frame(name="mean")
            .unstack("unit")
            .pint.quantify()
        )

        combined_df = pd.concat([mean_df, std_df], axis=1)
        groupby_levels = self.columns.names + no_x_levels

        for name, group in combined_df.groupby(groupby_levels):
            group = group.droplevel(groupby_levels)

            label = [item for pair in zip(groupby_levels, name) for item in pair]

            label = " : ".join(map(str, label))

            ax.plot(group.index, group['mean'], linewidth=1, linestyle='--', **kwargs)
            ax.fill_between(
                x=group.index,
                y1=group['mean'] + group['std'],
                y2=group['mean'] - group['std'],
                alpha=alpha,
                edgecolor="black",
                label=label
            )

        ax.legend()

    def _format_legend(self, ax: plt.Axes) -> None:
        leg = ax.get_legend()  # Get the existing legend from the axes
        for text in leg.get_texts():
            original_label = text.get_text()
            new_label = original_label.replace(')', '').replace('(', '').replace(', ', ' | ')
            text.set_text(new_label)

    def _plot_without_std(self, ax: plt.Axes, x: str, **kwargs) -> None:
        """
        Plots the data without standard deviation shading.
        Handles both real and imaginary parts for complex data.

        Parameters
        ----------
        ax : plt.Axes
            The matplotlib Axes on which to plot.
        x : str
            The index level to use for the x-axis.
        y : Optional[str], optional
            The column to use for the y-axis. If None, the first available column is used.
        show : bool, optional
            Whether to display the plot.
        log_scale_x : bool, optional
            If True, sets a logarithmic scale for the x-axis.
        log_scale_y : bool, optional
            If True, sets a logarithmic scale for the y-axis.
        **kwargs : dict
            Additional keyword arguments for the plot.

        Returns
        -------
        None
        """
        df = self.copy()

        groupby_levels = df._get_complementary_axis(x)

        df = df.unstack(groupby_levels)

        if isinstance(df.index, pd.MultiIndex) and df.index.nlevels == 1:
            df.index = df.index.get_level_values(0)

        super(PyMieSimDataFrame, df).plot(ax=ax, **kwargs)

        legend = ax.legend()

        legend.set_title(' | '.join(df.columns.names))

    def add_weight(self, weight_index: str, weight: np.ndarray) -> "PyMieSimDataFrame":

        self._validate_axis(axis=weight_index)

        stacking_index = self._get_complementary_axis(weight_index)

        return (
            self
            .unstack(stacking_index)
            .multiply(weight.squeeze(), axis='index')
            .stack(stacking_index)
        )

    def sum_over(self, axis: str) -> "PyMieSimDataFrame":

        self._validate_axis(axis=axis)

        stacking_index = [name for name in self.index.names if name != axis]

        return self.groupby(stacking_index).sum()#.to_frame().T

