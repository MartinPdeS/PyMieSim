"""
Sphere: Coupling vs numerical aperture
======================================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim import experiment
from PyMieSim import single
from PyMieSim.units import degree, watt, AU, RIU, meter
from PyOptik import Material

# %%
# Defining the source to be employed.
source = experiment.source.Gaussian(
    wavelength=50e-9 * meter,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = experiment.scatterer.Sphere(
    diameter=500e-9 * meter,
    material=Material.BK7,
    medium_index=1 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = experiment.detector.Photodiode(
    NA=np.linspace(0.1, 1.9, 1500) * AU,
    phi_offset=0 * degree,
    gamma_offset=0 * degree,
    polarization_filter=None,
    sampling=2000 * AU
)

# %%
# Defining the experiment setup
setup = experiment.Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

# %%
# Measuring the properties
data = setup.get(experiment.measure.coupling, export_as='dataframe')

# %%
# Plotting the results
figure = data.plot_data(x='NA')

single_source = single.Gaussian(
    wavelength=950e-9,
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)

single_scatterer = single.scatterer.Sphere(
    diameter=500e-9,
    material=Material.BK7,
    medium_index=1,
    source=single_source
)


print(single_scatterer.Qsca * 1e-3)