"""
Bohren-Huffman (figure~8.7)
===========================

"""

# %%
# Importing the dependencies: numpy, matplotlib, PyMieSim
import numpy
import matplotlib.pyplot as plt

from PyMieSim.tools.directories import validation_data_path

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure

theoretical = numpy.genfromtxt(f"{validation_data_path}/Figure87BH.csv", delimiter=',')

diameter = numpy.geomspace(10e-9, 6e-6, 800)
volume = numpy.pi * (diameter / 2)**2

source_set = Gaussian(
    wavelength=632.8e-9,
    polarization_value=[0, 90],
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

scatterer_set = Cylinder(
    diameter=diameter,
    index=1.55,
    n_medium=1,
    source_set=source_set
)

experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=None
)

data = experiment.get(measure.Csca)
data = data.y.values.squeeze() / volume * 1e-4 / 100

plt.figure(figsize=(8, 4))
plt.plot(diameter, data[0], 'b--', linewidth=3, label='PyMieSim')
plt.plot(diameter, data[1], 'r-', linewidth=3, label='PyMieSim')

plt.plot(diameter, theoretical[0], 'k--', linewidth=1, label='BH 8.8')
plt.plot(diameter, theoretical[1], 'k--', linewidth=1, label='BH 8.8')

plt.xlabel(r'Diameter [$\mu$m]')
plt.ylabel('Scattering cross section [cylinder]')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

assert numpy.all(numpy.isclose(data, theoretical, 1e-9)), 'Error: mismatch on BH_8.7 calculation occuring'

# -
