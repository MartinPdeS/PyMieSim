"""
PyMieSim vs PyMieScatt: 0
=========================

"""

# %%
# Importing the dependencies: numpy, matplotlib, PyMieSim
import numpy
import matplotlib.pyplot as plt

from PyMieSim.tools.directories import validation_data_path

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure

data_directory = validation_data_path.joinpath('PyMieScattQscaMedium.csv')
theoretical = numpy.genfromtxt(data_directory, delimiter=',')

diameter = numpy.geomspace(10e-9, 6e-6, 800)

source_set = Gaussian(
    wavelength=632.8e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)
scatterer_set = Sphere(
    diameter=diameter,
    index=1.4,
    n_medium=1.21,
    source_set=source_set
)


experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=None
)

data = experiment.get(measure.Qsca)
data = data.y.values.squeeze()

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
