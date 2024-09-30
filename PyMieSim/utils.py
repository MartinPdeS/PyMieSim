from MPSPlots.styles import mps
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional
import pint_pandas

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


from itertools import product
import numpy as np
from PyMieSim.units import Quantity
import pandas as pd
from PyMieSim.directories import validation_data_path

def get_pymiescatt_sphere_dataframe(
        wavelengths: Quantity,
        diameters: Quantity,
        indexes: Quantity,
        medium_indexes: Quantity,
        save_name: str = None) -> pd.DataFrame:
    """
    Generate a DataFrame with Mie scattering data for spheres, based on the provided
    wavelengths, diameters, and refractive indices. Optionally saves the DataFrame to CSV.

    Parameters:
    - wavelengths: Quantity array of wavelengths
    - diameters: Quantity array of sphere diameters
    - indexes: Quantity array of refractive indices
    - medium_indexes: Quantity array of medium refractive indices
    - save_name: Optional name for saving the DataFrame as a CSV file

    Returns:
    - pd.DataFrame: DataFrame containing Mie scattering properties.
    """
    import PyMieScatt as pms

    # Ensure the inputs are 1D arrays
    wavelengths = np.atleast_1d(wavelengths)
    diameters = np.atleast_1d(diameters)
    indexes = np.atleast_1d(indexes)
    medium_indexes = np.atleast_1d(medium_indexes)

    # Create a MultiIndex from the parameters
    indices = pd.MultiIndex.from_product([
            pint_pandas.PintArray(indexes, dtype=indexes.units),
            pint_pandas.PintArray(diameters, dtype=diameters.units),
            pint_pandas.PintArray(wavelengths, dtype=wavelengths.units),
            pint_pandas.PintArray(medium_indexes, dtype=medium_indexes.units),
        ],
        names=["index", "diameter", "wavelength", 'medium_index']
    )

    # Define the DataFrame structure
    dataframe = pd.DataFrame(
        columns=['Qext', 'Qsca', 'Qabs', 'g', 'Qpr', 'Qback', 'Qratio'],
        index=indices
    )


    # Calculate Mie scattering properties for each combination
    # print(indexes, diameters, wavelengths, medium_indexes)
    # dsa
    dataframe[:] = [
        pms.MieQ(m=index, wavelength=wavelength, diameter=diameter, nMedium=medium_index)
        for index, diameter, wavelength, medium_index in product(
            indexes.magnitude,
            diameters.magnitude,
            wavelengths.magnitude,
            medium_indexes.magnitude
        )
    ]

    # Assign the results to the DataFrame
    dataframe = dataframe.astype(float)

    # Save the DataFrame if a save_name is provided
    if save_name:
        save_path = validation_data_path / 'pymiescatt' / f"{save_name}.csv"
        dataframe.to_csv(save_path)

    return dataframe


def get_pymiescatt_coreshell_dataframe(
        wavelengths: Quantity,
        core_diameters: Quantity,
        shell_widths: Quantity,
        shell_indexes: Quantity,
        core_indexes: Quantity,
        medium_indexes: Quantity,
        save_name: str = None) -> pd.DataFrame:
    """
    Generate a DataFrame with Mie scattering data for spheres, based on the provided
    wavelengths, diameters, and refractive indices. Optionally saves the DataFrame to CSV.

    Parameters:
    - wavelengths: Quantity array of wavelengths
    - core_diameters: Quantity array of core diameters
    - shell_widths: Quantity array of shell diameters
    - core_indexes: Quantity array of core refractive indices
    - shell_indexes: Quantity array of shell refractive indices
    - medium_indexes: Quantity array of medium refractive indices
    - save_name: Optional name for saving the DataFrame as a CSV file

    Returns:
    - pd.DataFrame: DataFrame containing Mie scattering properties.
    """
    import PyMieScatt as pms

    # Ensure the inputs are 1D arrays
    wavelengths = np.atleast_1d(wavelengths)
    core_diameters = np.atleast_1d(core_diameters)
    shell_widths = np.atleast_1d(shell_widths)
    core_indexes = np.atleast_1d(core_indexes)
    shell_indexes = np.atleast_1d(shell_indexes)
    medium_indexes = np.atleast_1d(medium_indexes)

    # Create a MultiIndex from the parameters
    indices = pd.MultiIndex.from_product([
            pint_pandas.PintArray(core_indexes, dtype=core_indexes.units),
            pint_pandas.PintArray(shell_indexes, dtype=shell_indexes.units),
            pint_pandas.PintArray(core_diameters, dtype=core_diameters.units),
            pint_pandas.PintArray(shell_widths, dtype=shell_widths.units),
            pint_pandas.PintArray(wavelengths, dtype=wavelengths.units),
            pint_pandas.PintArray(medium_indexes, dtype=medium_indexes.units),
        ],
        names=["core_index", "shell_index", "core_diameter", "shell_width", "wavelength", 'medium_index']
    )

    # Define the DataFrame structure
    dataframe = pd.DataFrame(
        columns=['Qext', 'Qsca', 'Qabs', 'g', 'Qpr', 'Qback', 'Qratio'],
        index=indices
    )

    # Calculate Mie scattering properties for each combination
    dataframe[:] = [
        pms.MieQCoreShell(
            mCore=core_index,
            mShell=shell_index,
            wavelength=wavelength,
            dCore=core_diameter,
            dShell=core_diameter + shell_width,
            nMedium=medium_index
        )
        for core_index, shell_index, core_diameter, shell_width, wavelength, medium_index in product(
            core_indexes.to_base_units().magnitude,
            shell_indexes.to_base_units().magnitude,
            core_diameters.to_base_units().magnitude,
            shell_widths.to_base_units().magnitude,
            wavelengths.to_base_units().magnitude,
            medium_indexes.to_base_units().magnitude
        )
    ]

    dataframe = dataframe.astype(float)

    # Save the DataFrame if a save_name is provided
    if save_name:
        save_path = validation_data_path / 'pymiescatt' / f"{save_name}.csv"
        dataframe.to_csv(save_path)

    return dataframe