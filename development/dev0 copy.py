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


res = Material.polystyren.compute_refractive_index([780e-9, 900e-9])

print(res)

dsa

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=numpy.linspace(780, 900, 200) * nanometer,
    polarization=0 * degree,
    optical_power=10e-3 * watt,
    NA=[0.2] * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=[100, 120] * nanometer,
    property=Material.polystyren,
    medium_property=1.33 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number=['LP01'],
    NA=0.2 * AU,
    rotation=0 * degree,
    phi_offset=180.0 * degree,
    gamma_offset=0.0 * degree,
    sampling=600 * AU,
    mean_coupling=False,
    polarization_filter=None
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('coupling', drop_unique_level=True, scale_unit=True)

# %%
# Plotting the results
dataframe.plot_data(x='source:wavelength')
