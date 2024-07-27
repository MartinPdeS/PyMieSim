"""
Sphere: Coupling vs numerical aperture
======================================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim import experiment
from PyMieSim import single
from PyMieSim.experiment import measure
from PyOptik import UsualMaterial

# %%
# Defining the source to be employed.
source = experiment.source.Gaussian(
    wavelength=50e-9,
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = experiment.scatterer.Sphere(
    diameter=500e-9,
    material=UsualMaterial.BK7,
    medium_index=1,
    source=source
)

# %%
# Defining the detector to be employed.
detector = experiment.detector.Photodiode(
    NA=np.linspace(0.1, 1.9, 1500),
    phi_offset=0,
    gamma_offset=0,
    polarization_filter=[None],
    sampling=2000
)

# %%
# Defining the experiment setup
experiment = experiment.Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

# %%
# Measuring the properties
data = experiment.get(measure.coupling)

# %%
# Plotting the results
figure = data.plot(
    x=detector.NA,
)

single_source = single.Gaussian(
    wavelength=950e-9,
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)

single_scatterer = single.scatterer.Sphere(
    diameter=500e-9,
    material=UsualMaterial.BK7,
    medium_index=1,
    source=single_source
)


print(single_scatterer.Qsca * 1e-3)

_ = figure.show()
