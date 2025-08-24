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
from TypedUnit import ureg
import matplotlib.pyplot as plt

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup


permitivity = numpy.linspace(-10, 50, 400)

index = numpy.sqrt(permitivity.astype(complex)) * ureg.RIU

diameter = numpy.linspace(1, 200, 400) * ureg.nanometer

source = Gaussian(
    wavelength=400 * ureg.nanometer,
    polarization=90 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU
)


scatterer = Sphere(
    diameter=diameter,
    property=index,
    medium_property=1 * ureg.RIU,
    source=source
)

experiment = Setup(
    scatterer=scatterer,
    source=source
)

data = experiment.get('Qsca', add_units=False).squeeze().values.reshape([permitivity.size, diameter.size])

figure, ax = plt.subplots(1, 1)
ax.set(
    xlabel="Permittivity",
    ylabel=r'Diameter [$\mu$ m]',
    title="Scattering efficiency of a sphere"
)

image = ax.pcolormesh(permitivity, diameter, numpy.log(data))

plt.colorbar(mappable=image)

plt.show()
