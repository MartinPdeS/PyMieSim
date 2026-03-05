"""
CoreShell: Coupling vs Diameter
===============================

This example demonstrates how to compute and visualize the coupling efficiency as a function of core diameter for CoreShell scatterers using PyMieSim.
"""

import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import PhotodiodeSet
from PyMieSim.experiment.scatterer import CoreShellSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[1.2] * ureg.micrometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = CoreShellSet(
    core_diameter=numpy.geomspace(100, 600, 400) * ureg.nanometer,
    shell_thickness=[800] * ureg.nanometer,
    core_material=[Material.silver],
    shell_material=[Material.BK7],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1] * ureg.AU,
    phi_offset=[-180.0] * ureg.degree,
    gamma_offset=[0.0] * ureg.degree,
    sampling=[600],
    polarization_filter=[1] * ureg.degree,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="scatterer:core_diameter")
