from MPSPlots.styles import mps
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional

def generate_dataframe(experiment, is_complex: bool = False):
    """
    Generates a pandas DataFrame with a MultiIndex based on the mapping
    of 'source' and 'scatterer' from an experiment object.

    Parameters:
    experiment: The experiment object containing 'source' and 'scatterer' mappings.

    Returns:
    pd.DataFrame: A DataFrame with a MultiIndex created from the experiment mappings.
    """
    # Combine 'source' and 'scatterer' mappings for values and names
    iterables = {
        **experiment.source.mapping, **experiment.scatterer.mapping
    }

    if experiment.detector is not None:
        iterables.update(experiment.detector.mapping)

    if is_complex:
        iterables['type'] = ['real', 'imag']

    # Create a MultiIndex from the iterables
    row_index = pd.MultiIndex.from_product(iterables.values(), names=iterables.keys())

    # Return an empty DataFrame with the generated MultiIndex
    df = pd.DataFrame(index=row_index)

    df.columns.names = ['Data']

    return df


def plot_dataframe(
        dataframe: pd.DataFrame,
        x: str,
        alpha: float = 0.4,
        std: Optional[str] = None) -> None:
    """
    Plots a DataFrame using a specified MultiIndex level for the x-axis and
    optionally shades the mean ± standard deviation if provided.

    Parameters:
        dataframe (pd.DataFrame): The DataFrame to plot with a MultiIndex.
        x (str): The MultiIndex level to use for the x-axis.
        alpha (float): The transparency level for the shaded region (default is 0.4).
        std (Optional[str]): The MultiIndex level to use for calculating the standard deviation.  If None, no standard deviation will be plotted.

    """
    # Drop levels with only one unique value to avoid redundancy
    unique_levels = [
        l for l in dataframe.index.names if dataframe.index.get_level_values(l).nunique() == 1
    ]
    df_cleaned = dataframe.droplevel(unique_levels)

    # Unstack all levels except the one specified for the x-axis
    no_x_levels = [l for l in df_cleaned.index.names if l != x]

    df_unstacked = df_cleaned.unstack(level=no_x_levels)

    # Calculate mean values for plotting
    with plt.style.context(mps):
        _, ax = plt.subplots(1, 1)

    if std is not None:
        plot_with_std(df_unstacked, ax, x, std, alpha)

    else:
        plot_without_std(df_unstacked, ax, x)

    # Set axis labels with units (if available)
    xlabel_cleaned = x.replace('_', ' ').title()
    units = getattr(df_unstacked.index.values, 'units', '')

    # Construct xlabel
    xlabel = f"{xlabel_cleaned} [{units}]" if units else xlabel_cleaned
    ax.set_xlabel(xlabel)

    # Get the legend from the axes
    legend = ax.get_legend()

    # Capitalize the legend labels
    txt = legend.get_title().get_text().replace(',', ', ').replace('_', ' ').title()
    legend.get_title().set_text(txt)

    # Show the plot
    plt.show()

def plot_with_std(dataframe: pd.DataFrame, ax: plt.Axes, x: str, std: str, alpha: float) -> None:
    """
    Plots the mean with a shaded region representing mean ± standard deviation.

    Parameters:
        dataframe (pd.DataFrame): The DataFrame after unstacking.
        ax (plt.Axes): The matplotlib axis on which to plot.
        x (str): The index level for the x-axis.
        std (str): The index level for calculating standard deviation.
        alpha (float): Transparency for the shaded area.
    """
    no_std_levels = [l for l in dataframe.columns.names if l not in [x, std]]

    grouped = dataframe.T.groupby(level=no_std_levels)

    group_mean = grouped.mean().T.plot(ax=ax, linestyle='--', linewidth=1)

    for _, group in grouped:
        group_mean = group.mean()
        group_std = group.std()

        ax.fill_between(
            x=group_mean.index,
            y1=group_mean - group_std,
            y2=group_mean + group_std,
            alpha=alpha,
            edgecolor='black',
            label=f'Mean ± 1 Std ({std})'
        )

def plot_without_std(dataframe: pd.DataFrame, ax: plt.Axes, x: str) -> None:
    """
    Plots the data without standard deviation shading. Handles real and imaginary parts if the data is complex.

    Parameters:
        df_unstacked (pd.DataFrame): The DataFrame after unstacking.
        ax (plt.Axes): The matplotlib axis on which to plot.
        x (str): The index level for the x-axis.
    """
    no_x_levels = [level for level in dataframe.columns.names if level not in [None, x]]

    if len(no_x_levels) != 0:
        dataframe.groupby(no_x_levels, axis=1).plot(ax=ax)
    else:
        dataframe.plot(ax=ax)
