
import pandas as pd
import numpy as np
import pint_pandas
from PyMieSim.units import AU, pint
from itertools import product
from PyMieSim.units import Quantity
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

    dataframe[:] = [
        pms.MieQ(
            m=index,
            wavelength=wavelength,
            diameter=diameter,
            nMedium=medium_index
        )
        for index, diameter, wavelength, medium_index in product(
            indexes.magnitude,
            diameters.magnitude,
            wavelengths.magnitude,
            medium_indexes.magnitude
        )
    ]

    # Assign the results to the DataFrame
    for name, col in dataframe.items():
        dataframe[name] = pint.PintArray(col.values.astype(float), dtype=AU)

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

    for name, col in dataframe.items():
        dataframe[name] = pint.PintArray(col.values.astype(float), dtype=AU)

    # Save the DataFrame if a save_name is provided
    if save_name:
        save_path = validation_data_path / 'pymiescatt' / f"{save_name}.csv"
        dataframe.to_csv(save_path)

    return dataframe