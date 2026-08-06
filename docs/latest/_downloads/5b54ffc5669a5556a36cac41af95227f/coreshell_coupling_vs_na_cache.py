"""
CoreShell: Coupling vs Diameter
===============================

This example demonstrates how to compute and visualize the coupling efficiency as a function of core diameter for CoreShell scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim import (
    ureg,
    PhotodiodeSet,
    CoreShellSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    SellmeierMaterial,
    TabulatedMaterial,
    SellmeierMedium,
)


polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[0.5] * ureg.micrometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = CoreShellSet(
    core_diameter=[1000, 1250, 1500] * ureg.nanometer,
    shell_thickness=[800] * ureg.nanometer,
    core_material=[TabulatedMaterial("silver")],
    shell_material=[SellmeierMaterial("BK7")],
    medium=[SellmeierMedium("water")],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.3],
    cache_numerical_aperture=numpy.linspace(0.0, 0.2, 200),
    phi_offset=[-180.0] * ureg.degree,
    gamma_offset=[0.0] * ureg.degree,
    sampling=[1000],
    polarization_filter=[1] * ureg.degree,
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.get("coupling")

result.plot(x="detector:cache_NA")
