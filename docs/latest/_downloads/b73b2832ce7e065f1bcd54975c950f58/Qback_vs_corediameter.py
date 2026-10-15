"""
CoreShell: Qback vs Core Diameter
=======================================

This example demonstrates how to compute and visualize the backscattering efficiency (Qback)
as functions of core diameter for CoreShell scatterers using PyMieSim.
"""

import PyOptik
from PyMieSim import from_pyoptik

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim import (
    Measure,
    ureg,
    CoreShellSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)



polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[800, 900, 1000] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = CoreShellSet(
    core_diameter=numpy.geomspace(100, 900, 400) * ureg.nanometer,
    shell_thickness=[800] * ureg.nanometer,
    core_material=[from_pyoptik(PyOptik.load_material(name="main/Ag/Johnson"))],
    shell_material=[from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7"))],
    medium=[from_pyoptik(PyOptik.load_material(name="main/H2O/Daimon-19.0C"), medium=True)],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.get(Measure.QBACK)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.attrs["units"])

result.plot(x="scatterer:core_diameter")
