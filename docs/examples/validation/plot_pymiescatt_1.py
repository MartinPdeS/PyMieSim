"""
PyMieSim vs PyMieScatt: 1
=========================

"""

# %%
# Importing the dependencies: numpy, matplotlib, PyMieSim, PyMieScatt
import numpy
import matplotlib.pyplot as plt

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure

import PyMieScatt as ps

diameters = numpy.geomspace(10e-9, 6e-6, 800)
index = 1.4

source_set = Gaussian(
    wavelength=632.8e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1,
    NA=0.2
)

scatterer_set = Sphere(
    diameter=diameters,
    index=index,
    n_medium=1.,
    source_set=source_set
)

PyMieScatt_data = []
for diameter in diameters:
    efficiencies = ps.MieQ(
        m=index,
        wavelength=source_set.wavelength.values,
        diameter=diameter,
    )

    PyMieScatt_data.append(efficiencies[1])

PyMieScatt_data = numpy.asarray(PyMieScatt_data).squeeze()

experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=None
)

data = experiment.get(measure.Qsca)
PyMieSim_data = data.y.values.squeeze()

plt.figure(figsize=(8, 4))
plt.plot(
    diameters,
    PyMieSim_data,
    'C1-',
    linewidth=3,
    label='PyMieSim'
)

plt.plot(
    diameters,
    PyMieScatt_data,
    'k--',
    linewidth=1,
    label='PyMieScatt'
)

plt.xlabel(r'diameter [$\mu$m]')
plt.ylabel('Scattering efficiency [Sphere]')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# -
