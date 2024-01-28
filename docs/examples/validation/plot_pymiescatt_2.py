"""
PyMieSim vs PyMieScatt: 2
=========================

"""

# %%
# Importing the dependencies: numpy, matplotlib, PyMieSim, PyMieScatt
import numpy
import matplotlib.pyplot as plt
import PyMieScatt as ps

from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure

core_diameters = numpy.geomspace(10e-9, 500e-9, 400)
shell_width = 600e-9
core_index = 1.5
shell_index = 1.4

source_set = Gaussian(
    wavelength=600e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1,
    NA=0.2
)

scatterer_set = CoreShell(
    core_diameter=core_diameters,
    shell_width=shell_width,
    core_index=core_index,
    shell_index=shell_index,
    n_medium=1.0,
    source_set=source_set
)

PyMieScatt_data = []
for core_diameter in core_diameters:
    efficiencies = ps.MieQCoreShell(
        mCore=core_index,
        mShell=shell_index,
        wavelength=source_set.wavelength.values,
        dCore=core_diameter,
        dShell=core_diameter + shell_width
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
plt.plot(core_diameters, PyMieSim_data, 'C1-', linewidth=3, label='PyMieSim')

plt.plot(core_diameters, PyMieScatt_data, 'k--', linewidth=1, label='PyMieScatt')

plt.xlabel(r'diameter [$\mu$m]')
plt.ylabel('Scattering efficiency [CoreShell]')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# -
