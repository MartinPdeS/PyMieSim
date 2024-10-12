"""
Samples Properties
==================
"""

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.units import degree, nanometer, AU, RIU, watt

source = Gaussian(
    wavelength=1000 * nanometer,
    polarization=0 * degree,
    optical_power=1 * watt,
    NA=0.3 * AU
)

source.plot()

scatterer = Sphere(
    diameter=800 * nanometer,
    source=source,
    property=1.4 * RIU,
    medium_property=1.0 * RIU
)


farfield = scatterer.get_far_field(sampling=300)

farfield.plot()

scatterer.print_properties()

# Property        Value
# --------------  ----------------------
# size_parameter  2.513274122871834
# area            502654.8245743669 nm²
# g               0.6960064419371799
# Qsca            1.7376473038362086
# Qext            1.7376473038362086
# Qabs            0.0
# Qback           0.3550202921419265
# Qratio          0.20431090438096802
# Qpr             0.5282335865514354
# Csca            873436.800681911 nm²
# Cext            873436.800681911 nm²
# Cabs            0.0 m²
# Cback           178452.66266694054 nm²
# Cratio          102697.86180024572 nm²
# Cpr             265519.1607823004 nm²
