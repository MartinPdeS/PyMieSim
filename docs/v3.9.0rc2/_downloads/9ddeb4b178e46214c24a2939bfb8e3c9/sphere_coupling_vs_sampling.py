"""
Sphere: Coupling vs sampling
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from TypedUnit import ureg

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=400 * ureg.nanometer,
    polarization=90 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=5000 * ureg.nanometer,
    property=Material.BK7,
    medium_property=1 * ureg.RIU,
    source=source,
)

# %%
# Defining the detector to be employed.
detector = Photodiode(
    NA=[0.2] * ureg.AU,
    phi_offset=numpy.linspace(-20, 20, 400) * ureg.degree,
    gamma_offset=0 * ureg.degree,
    sampling=[20, 40, 80, 160] * ureg.AU,
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
