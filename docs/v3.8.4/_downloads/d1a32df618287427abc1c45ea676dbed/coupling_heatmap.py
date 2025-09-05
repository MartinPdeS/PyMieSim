"""
Coupling heatmap of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.
"""

import numpy
from TypedUnit import ureg
import matplotlib.pyplot as plt

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup


index = numpy.linspace(1.3, 2.1, 100) * ureg.RIU

diameter = numpy.linspace(1, 2000, 100) * ureg.nanometer

source = Gaussian(
    wavelength=400 * ureg.nanometer,
    polarization=90 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)


scatterer = Sphere(
    diameter=diameter, property=index, medium_property=1 * ureg.RIU, source=source
)


detector = Photodiode(
    polarization_filter=0 * ureg.degree,
    NA=0.3 * ureg.RIU,
    phi_offset=0 * ureg.degree,
    gamma_offset=0 * ureg.degree,
    sampling=400 * ureg.AU,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling", add_units=False)

values = dataframe.values.reshape([diameter.size, index.size])

figure, ax = plt.subplots(1, 1)

image = ax.pcolormesh(index, diameter * 1e6, values, shading="auto")
ax.set(
    xlabel="Permittivity", ylabel=r"Diameter [$\mu$m]", title="Coupling power [Watt]"
)
plt.colorbar(mappable=image)

plt.show()
