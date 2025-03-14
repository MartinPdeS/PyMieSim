"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Sphere, CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment.detector import Photodiode

from PyMieSim.experiment import Setup
from PyMieSim import units
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyOptik import Material
Material.print_available()



# %%
# Defining the source to be employed.
source_1 = Gaussian(
    wavelength=[405, 810] * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

scatterer_1 = Sphere(
    diameter=np.linspace(10, 1000, 150) * nanometer,
    medium_property=[1.33] * RIU,
    property=Material.polystyren,
    source=source_1
)

detector_1 = Photodiode(
    NA=[0.9 / 1.33, 1.3 / 1.33, 0.2, 0.1] * units.AU,
    gamma_offset=[0] * units.degree,
    phi_offset=90 * units.degree,

)

# %%
# Defining the experiment setup
experiment_1 = Setup(scatterer=scatterer_1, source=source_1, detector=detector_1)

dataframe_1 = experiment_1.get('coupling', scale_unit=True, drop_unique_level=True)

# dataframe_1._plot(x='scatterer:diameter', std='detector:NA', show=True)
dataframe_1._plot(x='scatterer:diameter', show=True)









# plt.show()