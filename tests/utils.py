import pandas as pd
import numpy as np

from itertools import product
from PyMieSim.units import Quantity, nanometer, RIU
from PyMieSim.directories import validation_data_path


def get_pymiescatt_sphere_dataframe(
        wavelength: Quantity,
        diameter: Quantity,
        index: Quantity,
        medium_index: Quantity,
        save_name: str = None) -> pd.DataFrame:
    """
    Generate a DataFrame with Mie scattering data for spheres, based on the provided
    wavelengths, diameters, and refractive indices. Optionally saves the DataFrame to CSV.

    Parameters
    ----------
    wavelength : Quantity
        array of wavelengths
    diameter : Quantity
        array of sphere diameter
    index : Quantity
        array of refractive indices
    medium_index : Quantity
        array of medium refractive indices
    save_name : Optional
        name for saving the DataFrame as a CSV file

    Returns
    -------
    pd.DataFrame
        DataFrame containing Mie scattering properties.
    """
    import PyMieScatt as pms

    # Ensure the inputs are 1D arrays
    wavelength = np.atleast_1d(wavelength).to_base_units().magnitude
    diameter = np.atleast_1d(diameter).to_base_units().magnitude
    index = np.atleast_1d(index).to_base_units().magnitude
    medium_index = np.atleast_1d(medium_index).to_base_units().magnitude

    # Create a MultiIndex from the parameters
    indices = pd.MultiIndex.from_product(
        [index, diameter, wavelength, medium_index],
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
            m=_index,
            wavelength=wavelength,
            diameter=diameter,
            nMedium=medium_index
        )
        for _index, diameter, wavelength, medium_index in product(
            index,
            diameter,
            wavelength,
            medium_index
        )
    ]

    # Assign the results to the DataFrame
    for name, col in dataframe.items():
        dataframe[name] = col.values

    # Save the DataFrame if a save_name is provided
    if save_name:
        save_path = validation_data_path / 'pymiescatt' / f"{save_name}.csv"
        dataframe.to_csv(save_path)

        print(f'Saving data: {save_path}')

    return dataframe


def get_pymiescatt_coreshell_dataframe(
        wavelength: Quantity,
        core_diameter: Quantity,
        shell_width: Quantity,
        shell_index: Quantity,
        core_index: Quantity,
        medium_index: Quantity,
        save_name: str = None) -> pd.DataFrame:
    """
    Generate a DataFrame with Mie scattering data for spheres, based on the provided
    wavelengths, diameters, and refractive indices. Optionally saves the DataFrame to CSV.

    Parameters
    ----------
    wavelength : Quantity
        array of wavelength
    core_diameter : Quantity
        array of core diameters
    shell_width : Quantity
        array of shell diameters
    core_index : Quantity
        array of core refractive indices
    shell_index : Quantity
        array of shell refractive indices
    medium_index : Quantity
        array of medium refractive indices
    save_name : Optional
        name for saving the DataFrame as a CSV file

    Returns
    -------
    pd.DataFrame
        DataFrame containing Mie scattering properties.
    """
    import PyMieScatt as pms

    # Ensure the inputs are 1D arrays
    wavelength = np.atleast_1d(wavelength).magnitude
    core_diameter = np.atleast_1d(core_diameter).magnitude
    shell_width = np.atleast_1d(shell_width).magnitude
    core_index = np.atleast_1d(core_index).magnitude
    shell_index = np.atleast_1d(shell_index).magnitude
    medium_index = np.atleast_1d(medium_index).magnitude

    # Create a MultiIndex from the parameters
    indices = pd.MultiIndex.from_product(
        [core_index, shell_index, core_diameter, shell_width, wavelength, medium_index],
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
            core_index,
            shell_index,
            core_diameter,
            shell_width,
            wavelength,
            medium_index
        )
    ]

    for name, col in dataframe.items():
        dataframe[name] = col.values

    # Save the DataFrame if a save_name is provided
    if save_name:
        save_path = validation_data_path / 'pymiescatt' / f"{save_name}.csv"
        dataframe.to_csv(save_path)

        print(f'Saving data: {save_path}')

    return dataframe




# Define sphere parameters
sphere_params = [
    dict(
        wavelength=632.8 * nanometer,
        medium_index=1.21 * RIU,
        index=1.4 * RIU,
        diameter=np.geomspace(10, 1_000, 50) * nanometer,
        save_name='example_sphere_0'
    ),
    dict(
        wavelength=632.8 * nanometer,
        medium_index=1.2 * RIU,
        index=(1.4 + 0.2j) * RIU,
        diameter=np.geomspace(10, 6_000, 800) * nanometer,
        save_name='example_sphere_1'
    )
]

# Define core-shell parameters
coreshell_params = [
    dict(
        wavelength=600 * nanometer,
        medium_index=1.0 * RIU,
        core_index=1.5 * RIU,
        shell_index=1.4 * RIU,
        shell_width=600 * nanometer,
        core_diameter=np.geomspace(10, 500, 40) * nanometer,
        save_name='example_coreshell_0'
    ),
    dict(
        wavelength=600 * nanometer,
        medium_index=1.0 * RIU,
        core_index=1.5 * RIU,
        shell_index=1.4 * RIU,
        shell_width=1200 * nanometer,
        core_diameter=np.geomspace(10, 500, 400) * nanometer,
        save_name='example_coreshell_1'
    )
]

# Compute dataframes for spheres
for params in sphere_params:
    get_pymiescatt_sphere_dataframe(**params)

# Compute dataframes for core-shell
for params in coreshell_params:
    get_pymiescatt_coreshell_dataframe(**params)
