"""
Sphere: Coupling vs wavelength
==============================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup



from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU


source = Gaussian(
    wavelength=np.linspace(950, 1050, 200) * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

scatterer = Sphere(
    diameter=100 * nanometer,
    property=Material.BK7,
    medium_property=Material.water,
    source=source
)

detector = CoherentMode(
    mode_number=["LP01"],
    NA=[0.05, 0.01] * AU,
    phi_offset=180 * degree,
    gamma_offset=0 * degree,
    polarization_filter=None,
    rotation=0 * degree,
    sampling=[300, 400, 500] * AU
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('coupling', scale_unit=True)

# %%
# Plotting the results
dataframe.plot_data(x="source:wavelength")
