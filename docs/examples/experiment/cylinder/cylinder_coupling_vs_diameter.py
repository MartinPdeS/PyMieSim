"""
Cylinder: Coupling vs Diameter
==============================

This example demonstrates how to compute and visualize the coupling efficiency as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

source = Gaussian(
    wavelength=[100, 1200] * ureg.nanometer,  # 1200 nm
    polarization=90 * ureg.degree,  # Polarization angle in ureg.degrees
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU,  # Numerical Aperture
)

scatterer = Sphere(
    diameter=np.linspace(100, 300, 200)
    * ureg.nanometer,  # Diameters ranging from 100 nm to 3000 nm
    property=[1.4] * ureg.RIU,  # Material of the cylinder
    medium_property=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source,
)

detector = Photodiode(
    NA=[0.1] * ureg.AU,  # Numerical Apertures for the detector
    phi_offset=[-180.0] * ureg.degree,  # Phi offset in ureg.degrees
    gamma_offset=[0.0] * ureg.degree,  # Gamma offset in ureg.degrees
    sampling=600 * ureg.AU,  # Number of sampling points
    polarization_filter=None,  # No polarization filter
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="scatterer:diameter")
