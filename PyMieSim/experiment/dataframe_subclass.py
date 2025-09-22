from typing import Self

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from MPSPlots import helper


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
            raise ValueError(
                f"x parameter '{axis}' is not in the DataFrame index. Available index levels: {available}"
            )

    def _get_complementary_axis(self, *axis) -> tuple:
        return [level for level in self.index.names if level not in axis]

    def plot(self, *args, **kwargs):
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
        figure_size : tuple
            Figure size.
        show : bool
            Plot the figure.
        save_as : str
            if set, save the figure.
        tight_layout : bool
            If set, render the plot in tight layout.
        show : bool, optional
            If True, displays the plot.
        **kwargs : dict
            Additional keyword arguments passed to the underlying plotting functions.

        Returns
        -------
        plt.Axes
            The matplotlib Axes with the plot.
        """
        if kwargs.get("std", False):
            return self.plot_with_std(*args, **kwargs)
        else:
            return self.plot_standard(*args, **kwargs)

    @helper.post_mpl_plot
    def plot_with_std(self, x: str, std: str, alpha: float = 0.5, **kwargs) -> None:
        """
        Plot the mean with standard deviation shading.
        Expects the DataFrame to have a MultiIndex and a 'pint' attribute.

        Parameters
        ----------
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
        figure, ax = plt.subplots(1, 1)

        self._validate_axis(axis=x)
        self._validate_axis(axis=std)

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

            ax.plot(group.index, group["mean"], linewidth=1, linestyle="--", **kwargs)
            ax.fill_between(
                x=group.index,
                y1=group["mean"] + group["std"],
                y2=group["mean"] - group["std"],
                alpha=alpha,
                edgecolor="black",
                label=label,
            )

        ax.legend()

        return figure

    def _format_legend(self, ax: plt.Axes) -> None:
        """
        Formats the legend of a Matplotlib Axes by replacing parentheses and commas
        in the legend labels with vertical bars for better readability.

        Parameters
        ----------
        ax : plt.Axes
            The Matplotlib Axes object whose legend will be formatted.
        """
        leg = ax.get_legend()  # Get the existing legend from the axes
        for text in leg.get_texts():
            original_label = text.get_text()
            new_label = (
                original_label.replace(")", "").replace("(", "").replace(", ", " | ")
            )
            text.set_text(new_label)

    @helper.post_mpl_plot
    def plot_standard(self, x: str, **kwargs) -> None:
        """
        Generate a line plot of the data without standard deviation shading.

        If the data is complex, both real and imaginary parts are plotted.

        Parameters
        ----------
        x : str
            Name of the index level to use for the x-axis.
        y : str, optional
            Name of the column to plot on the y-axis. If None, the first available column is used.
        show : bool, default=True
            If True, display the plot immediately.
        log_scale_x : bool, default=False
            If True, set the x-axis to logarithmic scale.
        log_scale_y : bool, default=False
            If True, set the y-axis to logarithmic scale.
        **kwargs
            Additional keyword arguments passed to the underlying Matplotlib plot call.

        Returns
        -------
        None
        """
        figure, ax = plt.subplots(1, 1)

        self._validate_axis(axis=x)

        df = self.copy()

        groupby_levels = df._get_complementary_axis(x)

        df = df.unstack(groupby_levels)

        if isinstance(df.index, pd.MultiIndex) and df.index.nlevels == 1:
            df.index = df.index.get_level_values(0)

        super(PyMieSimDataFrame, df).plot(ax=ax, **kwargs)

        legend = ax.legend()

        legend.set_title(" | ".join(df.columns.names))

        return figure

    def add_weight(self, weight_index: str, weight: np.ndarray) -> Self:
        """
        Multiply the DataFrame by a weight array along a specified MultiIndex level.

        Parameters
        ----------
        weight_index : str
            The MultiIndex level to apply the weights to.
        weight : np.ndarray
            The weight array to multiply with the DataFrame.

        Returns
        -------
        Self
            A new DataFrame with the weights applied.
        """
        self._validate_axis(axis=weight_index)

        stacking_index = self._get_complementary_axis(weight_index)

        return (
            self.unstack(stacking_index)
            .multiply(weight.squeeze(), axis="index")
            .stack(stacking_index)
        )

    def sum_over(self, axis: str) -> Self:
        """
        Sum the DataFrame over a specified MultiIndex level.
        Parameters
        ----------
        axis : str
            The MultiIndex level to sum over.

        Returns
        -------
        Self
            A new DataFrame summed over the specified axis.
        """
        self._validate_axis(axis=axis)

        stacking_index = [name for name in self.index.names if name != axis]

        return self.groupby(stacking_index).sum()  # .to_frame().T
