"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Sphere, CoreShell
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU, volt, meter


source = Gaussian(
    wavelength=[600, 700, 800] * nanometer,
    polarization=30 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)


# source = PlaneWave(
#     wavelength=[500, 600] * nanometer,
#     polarization=30 * degree,
#     amplitude=1e-3 * volt / meter,
# )

scatterer = Sphere(
    diameter=np.geomspace(6.36, 1000, 150) * nanometer,
    medium_property=1 * RIU,
    property=1.6 * RIU,
    source=source
)

# scatterer = CoreShell(
#     core_diameter=np.geomspace(6.36, 1000, 150) * nanometer,
#     shell_width=100 * nanometer,
#     medium_property=1 * RIU,
#     core_property=1.6 * RIU,
#     shell_property=1.6 * RIU,
#     source=source
# )


detector = Photodiode(
    NA=0.1 * AU,
    phi_offset=0 * degree,
    gamma_offset=0 * degree,
    sampling=100 * AU,
    polarization_filter=None,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)


dataframe = experiment.get('coupling', scale_unit=True, drop_unique_level=True)


dataframe.plot_data(x='diameter')