"""
Sphere: Coupling vs diameter
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200 * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=[0.1] * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=numpy.linspace(100, 10000, 600) * nanometer,
    property=Material.BK7,
    medium_property=1.0 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number=['LP01', 'LP11', 'LP02'],
    NA=0.2 * AU,
    rotation=0 * degree,
    phi_offset=0.0 * degree,
    gamma_offset=0.0 * degree,
    sampling=600 * AU,
    mean_coupling=True,
    polarization_filter=None
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('a1', drop_unique_level=True, scale_unit=True)

# %%
# Plotting the results
dataframe.plot(x='scatterer:diameter')
