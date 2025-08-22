"""
Sphere: Qsca vs wavelength
==========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import MaterialBank
# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=np.linspace(400, 1000, 50) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=[200] * ureg.nanometer,
    property=MaterialBank.BK7,
    medium_property=1 * ureg.RIU,
    source=source
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get('Qsca', 'Qpr', scale_unit=True)

# %%
# Plotting the results
dataframe.plot(x="source:wavelength")
