"""
PyMieSim vs PyMieScatt: 0
=========================

"""

# %%
# Importing the dependencies: numpy, matplotlib, PyMieSim
import numpy
import matplotlib.pyplot as plt

from PyMieSim.tools.directories import validation_data_path
from PyMieSim.experiment import SourceSet, SphereSet, Setup
from PyMieSim import measure

data_directory = validation_data_path.joinpath('PyMieScattQscaMedium.csv')
theoretical = numpy.genfromtxt(data_directory, delimiter=',')

diameter = numpy.geomspace(10e-9, 6e-6, 800)

scatterer_set = SphereSet(
    diameter=diameter,
    index=1.4,
    n_medium=1.21
)

source_set = SourceSet(
    wavelength=632.8e-9,
    linear_polarization=0,
    amplitude=1
)

experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=None
)

data = experiment.Get(measure.Qsca)
data = data.array.squeeze()

plt.figure(figsize=(8, 4))
plt.plot(
    diameter,
    data,
    'C1-',
    linewidth=3,
    label='PyMieSim'
)

plt.plot(
    diameter,
    theoretical,
    'k--',
    linewidth=1,
    label='PyMieScatt'
)

plt.xlabel(r'diameter [$\mu$m]')
plt.ylabel('Scattering efficiency [Sphere + n_medium]')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# -
