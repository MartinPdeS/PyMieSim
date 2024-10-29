"""
Sphere: Coupling vs wavelength
==============================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU
import numpy


size = 4
ones = numpy.ones(size)

wavelength = numpy.linspace(100, 400, size) * nanometer
medium_property = 1.0 * RIU
property = 1.5 * RIU
diameter = 100 * nanometer
polarization = 0 * degree
optical_power = 200 * watt
NA = 0.2 * AU

source_sequential = Gaussian(
    wavelength=wavelength,
    polarization=ones * polarization,
    optical_power=ones * optical_power,
    NA=ones * NA
)

scatterer_sequential = Sphere(
    diameter=ones * diameter,
    property=ones * property,
    medium_property=ones * medium_property,
    source=source_sequential
)

experiment_sequential = Setup(scatterer=scatterer_sequential, source=source_sequential)

values_sequential = experiment_sequential.get_sequential('Qsca')

print(values_sequential)

source = Gaussian(
    wavelength=wavelength,
    polarization=polarization,
    optical_power=optical_power,
    NA=NA
)

scatterer = Sphere(
    diameter=diameter,
    property=property,
    medium_property=medium_property,
    source=source
)

experiment = Setup(scatterer=scatterer, source=source)

values = experiment.get('Qsca', add_units=False).values

print(values.squeeze())
