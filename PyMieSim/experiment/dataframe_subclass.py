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
        # Validate that 'x' is a valid index level.
        if x not in self.index.names:
            available = ", ".join(self.index.names)
            raise ValueError(f"x parameter '{x}' is not in the DataFrame index. Available index levels: {available}")

        # Validate that 'x' is a valid index level.
        if std is not None and std not in self.index.names:
            available = ", ".join(self.index.names)
            raise ValueError(f"x parameter '{std}' is not in the DataFrame index. Available index levels: {available}")

        if ax is None:
            with plt.style.context(mps):
                _, ax = plt.subplots()
                ax.set(xscale=xscale, yscale=yscale)

        if std is not None:
            self._plot_with_std(ax=ax, x=x, std=std, alpha=alpha, **kwargs)
        else:
            self._plot_without_std(ax=ax, x=x, **kwargs)

        # self._format_legend(ax)

        if show:
            plt.show()

        return ax

    def _plot_with_std(
        self,
        ax: plt.Axes,
        x: str,
        std: str,
        alpha: float = 0.5,
        show: bool = True,
        **kwargs) -> None:
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
        no_std_levels = [level for level in self.index.names if level != std]
        no_x_levels = [level for level in self.index.names if level not in [x, std]]

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

        title = leg.get_title().get_text()
        leg.set_title(title.replace(',', ' | '))

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
        # Stack levels if necessary for proper grouping
        df_to_plot = self.copy()
        if 'type' in self.columns.names:
            df_to_plot = df_to_plot.stack('type', future_stack=True)
        df_to_plot = df_to_plot.stack('data', future_stack=True)

        groupby_levels = [level for level in df_to_plot.index.names if level != x]

        if groupby_levels:
            for name, group in df_to_plot.groupby(groupby_levels, dropna=False):
                group = group.droplevel(groupby_levels)

                label = [item for pair in zip(groupby_levels, name) for item in pair]

                label = " : ".join(map(str, label))

                values = group.squeeze().values
                if hasattr(values, 'quantity'):
                    ax.plot(group.index, values.quantity, label=label, **kwargs)
                else:
                    ax.plot(group.index, values, label=label, **kwargs)

        else:
            # If no groupby, determine the column to plot
            col = df_to_plot.columns[0]
            ax.plot(df_to_plot.index, df_to_plot[col], label=str(col), **kwargs)

        ax.legend()
        # groupby_levels = [level for level in self.index.names if level != x]

        # df = self.unstack(groupby_levels)
        # super(PyMieSimDataFrame, self).plot(ax=ax, **kwargs)

        # ax.set_xlabel(f"{x} [{df.index.values.units}]")

