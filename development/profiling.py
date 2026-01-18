"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.single.scatterer import Sphere, CoreShell
from PyMieSim.single.source import Gaussian


from PyMieSim import units
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyOptik import Material

Material.print_available()

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=[405] * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU,
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=100 * nanometer,
    medium_property=[1.33] * RIU,
    property=Material.polystyren,
    source=source,
)

# scatterer = CoreShell(
#     # core_diameter=np.linspace(10, 1000, 1500) * nanometer,
#     core_diameter=160 * nanometer,
#     shell_thickness=6 * units.nanometer,
#     medium_property=1.33 * RIU,
#     core_property=1.38 * units.RIU,
#     shell_property=1.48 * units.RIU,
#     source=source
# )

scatterer.get_spf().plot()
