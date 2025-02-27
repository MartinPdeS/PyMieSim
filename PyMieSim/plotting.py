from MPSPlots.styles import mps
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional





class PyMieSimDataFrame(pd.DataFrame):
    """
    A subclass of pandas DataFrame with a custom plot method.
    """

    @property
    def _constructor(self):
        """
        Ensures that operations returning DataFrames return instances of CustomDataFrame.
        """
        return PyMieSimDataFrame

    def plot(self, x: str, alpha: float = 0.4, std: Optional[str] = None, ax: plt.Axes = None, **kwargs) -> None:
        """
        Plots a DataFrame using a specified MultiIndex level for the x-axis and
        optionally shades the mean ± standard deviation if provided.

        Parameters
        ----------
        dataframe : pd.DataFrame
            The DataFrame to plot with a MultiIndex.
        x : str
            The MultiIndex level to use for the x-axis.
        alpha : float
            The transparency level for the shaded region (default is 0.4).
        std : Optional[str]
            The MultiIndex level to use for calculating the standard deviation.  If None, no standard deviation will be plotted.

        """
        # Calculate mean values for plotting
        if not ax:
            with plt.style.context(mps):
                _, ax = plt.subplots(1, 1)

        if std is not None:
            return self.plot_with_std(ax=ax, x=x, std=std, alpha=alpha, **kwargs)

        else:
            return self.plot_without_std(ax=ax, x=x, **kwargs)

    def plot_with_std(self, ax, x: str, std: str, alpha: float = 0.5, show: bool = True, **kwargs) -> None:
        """
        Plot the mean with standard deviation shading for a given dataframe.

        Parameters
        ----------
        dataframe :  pd.DataFrame
            A pandas DataFrame containing the data with pint units.
        x : str
            The index level used as the x-axis for the plot.
        std : str
            The index level representing the standard deviation group.
        alpha : float, optional
            The transparency level for the shaded area representing standard deviation.
            Default is 0.5.
        """

        # Identify levels excluding the 'std' level
        no_std_levels = [level for level in self.index.names if level != std]
        no_x_levels = [level for level in self.index.names if level not in [x, std]]

        # Calculate the mean and standard deviation, preserving pint units
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

        # Concatenate mean and std into a single DataFrame
        combined_df = pd.concat([mean_df, std_df], axis=1)

        # Determine grouping by column names and other levels excluding x and std
        groupby_levels = self.columns.names + no_x_levels

        for name, group in combined_df.groupby(groupby_levels):
            group = group.droplevel(groupby_levels)
            name_str = " : ".join(map(str, name))  # Join group names into a label

            # Plot the mean line
            ax.plot(group.index, group['mean'], linewidth=1, linestyle='--', **kwargs)

            # Shade the area representing the standard deviation
            ax.fill_between(
                x=group.index,
                y1=group['mean'] + group['std'],
                y2=group['mean'] - group['std'],
                alpha=alpha,
                edgecolor="black",
                label=name_str
            )

        ax.legend(title=" : ".join(groupby_levels))

        if show:
            plt.show()


    def plot_without_std(
        self,
        ax: plt.Axes,
        x: str,
        y: str = None,
        show: bool = True,
        log_scale_x: bool = False,
        log_scale_y: bool = False,
        **kwargs) -> plt.Axes:
        """
        Plots data without standard deviation shading, handling both real and imaginary parts for complex data.

        Parameters
        ----------
        ax : plt.Axes
            The matplotlib axis on which to plot the data.
        x : str
            The column or index level to use as the x-axis.
        y : str, optional
            The column to use as the y-axis. If not specified, the first column is used.
        show : bool, default=True
            Whether to display the plot immediately.
        log_scale_x : bool, default=False
            Whether to apply logarithmic scaling to the x-axis.
        log_scale_y : bool, default=False
            Whether to apply logarithmic scaling to the y-axis.
        **kwargs : dict, optional
            Additional keyword arguments passed to `matplotlib.pyplot.plot`,
            such as line styles, colors, markers, etc.

        Returns
        -------
        plt.Axes
            The matplotlib axis with the plotted data.

        Notes
        -----
        - If the `dataframe` contains complex data, the function plots the real part by default.
        - The `kwargs` can be used to customize the appearance of the plot.
        - For multi-indexed DataFrames, groups are automatically handled and labeled.
        """
        if 'type' in self.columns.names:
            dataframe = self.stack('type', future_stack=True)
        else:
            dataframe = self

        dataframe = dataframe.stack('data', future_stack=True)

        groupby_levels = [level for level in dataframe.index.names if level not in [x]]

        if groupby_levels:
            for name, group in dataframe.groupby(groupby_levels, dropna=False):
                group = group.droplevel(groupby_levels)
                label = " : ".join(map(str, name))

                values = group.squeeze().values
                if hasattr(values, 'quantity'):
                    ax.plot(group.index, values.quantity, label=label)
                else:
                    ax.plot(group.index, values, label=label)

        else:
            name = dataframe.columns[0]
            label = " : ".join(map(str, name))
            ax.plot(
                dataframe.index,
                dataframe.squeeze(),
                label=label,
                **kwargs
            )

        # Apply log scaling if specified
        if log_scale_x:
            ax.set_xscale('log')
        if log_scale_y:
            ax.set_yscale('log')

        ax.legend(title=" : ".join(groupby_levels))

        if show:
            plt.show()

        return ax
