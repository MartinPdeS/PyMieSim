"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
import matplotlib.pyplot as plt
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyOptik import Material

from PyMieSim.experiment import Setup
from PyMieSim import units
from PyMieSim.units import degree, watt, AU, RIU

c = 2.998e8 * units.meter / units.second  # Speed of light in vacuum

start_frequency = 0.1 * units.gigahertz
end_frequency = 6 * units.gigahertz

frequencies = np.linspace(start_frequency, end_frequency, 250)

wavelengths = c / frequencies


source = Gaussian(
    wavelength=wavelengths,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU,
)

scatterer = Sphere(
    diameter=14e7 * units.nanometer,
    medium_property=[1] * RIU,
    property=Material.silver,
    source=source,
)


# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

data = experiment.get(
    "Qforward", "Qback", scale_unit=False, drop_unique_level=True, as_numpy=True
)

q_forward, q_back = data

plt.figure(figsize=(10, 6))
plt.plot(frequencies, q_forward, label="Qforward", color="blue")
plt.plot(frequencies, q_back, label="Qback", color="orange")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Scattering Efficiency")
plt.title("Scattering Efficiency vs Wavelength")
plt.legend()
plt.grid()
plt.yscale("log")
plt.tight_layout()
plt.show()
