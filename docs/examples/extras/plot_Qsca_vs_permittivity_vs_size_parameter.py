"""
Scattering efficiency of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.
"""

# %%
# Importing the package: PyMieSim
import numpy

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure
from PyMieSim.units import degree, watt, AU, nanometer, RIU

import matplotlib.pyplot as plt
from MPSPlots.render2D import SceneList


permitivity = numpy.linspace(-10, 50, 400)

index = numpy.sqrt(permitivity.astype(complex)) * RIU

diameter = numpy.linspace(1, 200, 400) * nanometer

source = Gaussian(
    wavelength=400 * nanometer,
    polarization=90 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)


scatterer = Sphere(
    diameter=diameter,
    index=index,
    medium_index=1 * RIU,
    source=source
)

experiment = Setup(
    scatterer=scatterer,
    source=source
)

data = experiment.get(measure.Qsca, export_as='numpy')

data = abs(data.squeeze())

figure, ax = plt.subplots(1, 1)
ax.set(
    xlabel="Permittivity",
    ylabel=r'Diameter [$\mu$ m]',
    title="Scattering efficiency of a sphere"
)

image = ax.pcolormesh(permitivity, diameter, numpy.log(data))

plt.colorbar(mappable=image)

plt.show()
