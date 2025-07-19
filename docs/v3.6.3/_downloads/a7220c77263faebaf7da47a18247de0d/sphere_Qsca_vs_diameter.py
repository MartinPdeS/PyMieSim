"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian

from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyOptik import Material
Material.print_available()
import PyMieSim
PyMieSim.debug_mode = True



# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=[405, 810] * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

scatterer = Sphere(
    diameter=np.linspace(10, 1000, 150) * nanometer,
    medium_property=[1.33] * RIU,
    property=Material.polystyren,
    source=source
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get('Qsca', scale_unit=True)

dataframe.plot(x='scatterer:diameter', show=True)









# plt.show()