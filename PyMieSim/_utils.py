import pandas as pd
import numpy as np

from itertools import product
from TypedUnit import ureg, AnyUnit
from PyMieSim.directories import validation_data_path


def get_pymiescatt_sphere_dataframe(
    wavelength: AnyUnit,
    diameter: AnyUnit,
    index: AnyUnit,
    medium_index: AnyUnit,
    save_name: str = None,
) -> pd.DataFrame:
    """
    Generate a DataFrame with Mie scattering data for spheres, based on the provided
    wavelengths, diameters, and refractive indices. Optionally saves the DataFrame to CSV.

    Parameters
    ----------
    wavelength : AnyUnit
        array of wavelengths
    diameter : AnyUnit
        array of sphere diameter
    index : AnyUnit
        array of refractive indices
    medium_index : AnyUnit
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
        names=["index", "diameter", "wavelength", "medium_index"],
    )

    # Define the DataFrame structure
    dataframe = pd.DataFrame(
        columns=["Qext", "Qsca", "Qabs", "g", "Qpr", "Qback", "Qratio"], index=indices
    )

    # Calculate Mie scattering properties for each combination
    dataframe[:] = [
        pms.MieQ(
            m=_index, wavelength=wavelength, diameter=diameter, nMedium=medium_index
        )
        for _index, diameter, wavelength, medium_index in product(
            index, diameter, wavelength, medium_index
        )
    ]

    # Assign the results to the DataFrame
    for name, col in dataframe.items():
        dataframe[name] = col.values

    # Save the DataFrame if a save_name is provided
    if save_name:
        save_path = validation_data_path / "pymiescatt" / f"{save_name}.csv"
        dataframe.to_csv(save_path)

        print(f"Saving data: {save_path}")

    return dataframe


def get_pymiescatt_coreshell_dataframe(
    wavelength: AnyUnit,
    core_diameter: AnyUnit,
    shell_width: AnyUnit,
    shell_index: AnyUnit,
    core_index: AnyUnit,
    medium_index: AnyUnit,
    save_name: str = None,
) -> pd.DataFrame:
    """
    Generate a DataFrame with Mie scattering data for spheres, based on the provided
    wavelengths, diameters, and refractive indices. Optionally saves the DataFrame to CSV.

    Parameters
    ----------
    wavelength : AnyUnit
        array of wavelength
    core_diameter : AnyUnit
        array of core diameters
    shell_width : AnyUnit
        array of shell diameters
    core_index : AnyUnit
        array of core refractive indices
    shell_index : AnyUnit
        array of shell refractive indices
    medium_index : AnyUnit
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
        names=[
            "core_index",
            "shell_index",
            "core_diameter",
            "shell_width",
            "wavelength",
            "medium_index",
        ],
    )

    # Define the DataFrame structure
    dataframe = pd.DataFrame(
        columns=["Qext", "Qsca", "Qabs", "g", "Qpr", "Qback", "Qratio"], index=indices
    )

    # Calculate Mie scattering properties for each combination
    dataframe[:] = [
        pms.MieQCoreShell(
            mCore=core_index,
            mShell=shell_index,
            wavelength=wavelength,
            dCore=core_diameter,
            dShell=core_diameter + shell_width,
            nMedium=medium_index,
        )
        for core_index, shell_index, core_diameter, shell_width, wavelength, medium_index in product(
            core_index,
            shell_index,
            core_diameter,
            shell_width,
            wavelength,
            medium_index,
        )
    ]

    for name, col in dataframe.items():
        dataframe[name] = col.values

    # Save the DataFrame if a save_name is provided
    if save_name:
        save_path = validation_data_path / "pymiescatt" / f"{save_name}.csv"
        dataframe.to_csv(save_path)

        print(f"Saving data: {save_path}")

    return dataframe


# Define sphere parameters
sphere_params = [
    dict(
        wavelength=632.8 * ureg.nanometer,
        medium_index=1.21 * ureg.RIU,
        index=1.4 * ureg.RIU,
        diameter=np.geomspace(10, 1_000, 50) * ureg.nanometer,
        save_name="example_sphere_0",
    ),
    dict(
        wavelength=632.8 * ureg.nanometer,
        medium_index=1.2 * ureg.RIU,
        index=(1.4 + 0.2j) * ureg.RIU,
        diameter=np.geomspace(10, 6_000, 800) * ureg.nanometer,
        save_name="example_sphere_1",
    ),
]

# Define core-shell parameters
coreshell_params = [
    dict(
        wavelength=600 * ureg.nanometer,
        medium_index=1.0 * ureg.RIU,
        core_index=1.5 * ureg.RIU,
        shell_index=1.4 * ureg.RIU,
        shell_width=600 * ureg.nanometer,
        core_diameter=np.geomspace(10, 500, 40) * ureg.nanometer,
        save_name="example_coreshell_0",
    ),
    dict(
        wavelength=600 * ureg.nanometer,
        medium_index=1.0 * ureg.RIU,
        core_index=1.5 * ureg.RIU,
        shell_index=1.4 * ureg.RIU,
        shell_width=1200 * ureg.nanometer,
        core_diameter=np.geomspace(10, 500, 400) * ureg.nanometer,
        save_name="example_coreshell_1",
    ),
]

# Compute dataframes for spheres
for params in sphere_params:
    get_pymiescatt_sphere_dataframe(**params)

# Compute dataframes for core-shell
for params in coreshell_params:
    get_pymiescatt_coreshell_dataframe(**params)
