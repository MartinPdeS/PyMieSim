"""
Sphere: Goniometer
==================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from TypedUnit import ureg

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.single.polarization import RightCircular
from PyOptik import Material

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200 * ureg.nanometer,
    polarization=RightCircular(),
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=2000 * ureg.nanometer,
    property=Material.BK7,
    medium_property=1 * ureg.RIU,
    source=source,
)

# %%
# Defining the detector to be employed.
detector = Photodiode(
    NA=[0.1] * ureg.AU,
    phi_offset=numpy.linspace(-180, 180, 200) * ureg.degree,
    cache_NA=0.05 * ureg.AU,
    gamma_offset=40 * ureg.degree,
    sampling=400 * ureg.AU,
    polarization_filter=[0, 30] * ureg.degree,
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get("coupling")

# %%
# Plotting the results
dataframe.plot(x="detector:phi_offset")
