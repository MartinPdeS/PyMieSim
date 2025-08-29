from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import plot_system
from math import sqrt

from PyMieSim.units import nanometer, degree, watt, AU, RIU


source = Gaussian(
    wavelength=1550 * nanometer,  # Wavelength of 1550 nm
    polarization=0 * degree,  # Linear polarization at 0 radians
    optical_power=1 * watt,  # Optical power in arbitrary units
    NA=0.3 * AU,  # Numerical Aperture
)


scatterer = Cylinder(
    diameter=7800 * nanometer,  # Diameter of 7.8 micrometers
    source=source,  # The Gaussian source defined above
    medium_property=1.0 * RIU,  # Refractive index of the surrounding medium
    property=sqrt(1.5) * RIU,  # Refractive index of the scatterer
)

detector = CoherentMode(
    mode_number="LP01",
    NA=0.2 * AU,  # Numerical Aperture
    gamma_offset=0 * degree,  # Gamma offset in degrees
    phi_offset=60 * degree,  # Phi offset in degrees
    polarization_filter=0 * degree,  # Polarization filter angle in degrees
)

spf = scatterer.get_spf()

# %%
# Visualize the system using the plot_system function
plot_system(scatterer, source, detector)
