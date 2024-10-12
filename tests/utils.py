from PyMieSim.utils import get_pymiescatt_sphere_dataframe, get_pymiescatt_coreshell_dataframe
import numpy as np
from PyMieSim.units import nanometer, RIU

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
