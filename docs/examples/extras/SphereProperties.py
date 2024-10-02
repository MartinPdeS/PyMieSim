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
    index=1.4 * RIU,
    medium_index=1.0 * RIU
)


farfield = scatterer.get_far_field(sampling=300)

farfield.plot()

scatterer.print_properties()

# Property              value
# --------------  -----------
# size_parameter  2.51327
# area            5.02655e-13
# index           1.4
# Qsca            1.73765
# Qext            1.73765
# Qabs            0
# Qback           0.35502
# Qratio          0.204311
# Qpr             0.528234
# Csca            8.73437e-13
# Cext            8.73437e-13
# Cabs            0
# Cback           1.78453e-13
# Cratio          0.204311
# Cpr             2.65519e-13
# g               0.696006
