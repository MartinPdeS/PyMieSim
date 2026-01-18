"""
CoreShell: Coupling vs Diameter
===============================

This example demonstrates how to compute and visualize the coupling efficiency as a function of core diameter for CoreShell scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

source = Gaussian(
    wavelength=1.2 * ureg.micrometer,  # 1200 nm
    polarization=90 * ureg.degree,  # Polarization angle in ureg.degrees
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU,  # Numerical Aperture
)

scatterer = CoreShell(
    core_diameter=[1500] * ureg.nanometer,  # Core diameters from 100 nm to 600 nm
    shell_thickness=800 * ureg.nanometer,  # Shell width of 800 nm
    core_property=Material.silver,  # Core material
    shell_property=Material.BK7,  # Shell material
    medium_property=1 * ureg.RIU,  # Surrounding medium's refractive index
    source=source,
)

detector = Photodiode(
    NA=[0.3] * ureg.AU,  # Numerical Apertures for the detector
    cache_NA=numpy.linspace(0.0, 0.2, 200) * ureg.AU,
    phi_offset=-180.0 * ureg.degree,  # Phi offset in ureg.degrees
    gamma_offset=0.0 * ureg.degree,  # Gamma offset in ureg.degrees
    sampling=4000 * ureg.AU,  # Number of sampling points
    polarization_filter=1 * ureg.degree,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:cache_NA")
