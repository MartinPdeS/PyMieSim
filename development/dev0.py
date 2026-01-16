import pytest
import numpy as np
from PyMieSim.units import ureg
from PyOptik import Material

from PyMieSim.single.detector import IntegratingSphere, Photodiode
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single import SystemPlotter

plotter = SystemPlotter()


source = Gaussian(
    wavelength=750 * ureg.nanometer,  # 750 nm
    polarization=30 * ureg.degree,  # Polarization in ureg.degrees
    optical_power=1 * ureg.watt,  # Power in ureg.watts
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

# Define the scatterer (sphere)
scatterer = Sphere(
    diameter=1 * ureg.nanometer,  # 1500 nm diameter
    source=source,
    property=1.8 * ureg.RIU,
    medium_property=1. * ureg.RIU
)

detector = IntegratingSphere(
    sampling=5000,
)


# 50.18021588708274 attowat

# plotter.plot(detector)

coupling = detector.get_coupling(scatterer=scatterer)
print("coupling", coupling)


# coupling 21.23395750875811 attowatt


# coupling 21.23395750875811 attowatt