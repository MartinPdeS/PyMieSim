"""
InfiniteCylinder: Goniometer
============================

This example demonstrates how to use a goniometer setup to measure and visualize the coupling efficiency as a function of angular displacement for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import PhotodiodeSet
from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = InfiniteCylinderSet(
    diameter=[2000] * ureg.nanometer,
    material=[Material.BK7],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

detector = PhotodiodeSet(
    numerical_aperture=[0.5, 0.3, 0.1, 0.05] * ureg.AU,
    phi_offset=np.linspace(-180, 180, 200) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[400],
    polarization_filter=None,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling", scale_unit=True)

dataframe.plot(x="detector:phi_offset")
