"""
Coupling heatmap of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.
"""

import numpy
from PyMieSim.units import ureg
import matplotlib.pyplot as plt

from PyMieSim.experiment.detector import PhotodiodeSet
from PyMieSim.experiment.scatterer import SphereSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup

from PyOptik import Material
import PyMieSim
PyMieSim.debug_mode = True

polarization_state = PolarizationSet(
    angles=[90] * ureg.degree,  # Linear polarization at 90 degrees
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_state,  # Linear polarization at 90 degrees
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = SphereSet(
    diameter=numpy.linspace(1, 2000, 100) * ureg.nanometer,
    refractive_index=numpy.linspace(1.3, 2.1, 100) * ureg.RIU,
    medium_material=[Material.water],
    source=source
)

detector = PhotodiodeSet(
    polarization_filter=[0] * ureg.degree,
    numerical_aperture=[0.3] * ureg.RIU,
    phi_offset=[0] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[400]
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

values = experiment.get("coupling", add_units=False, as_numpy=True)


figure, ax = plt.subplots(1, 1)

image = ax.pcolormesh(
    scatterer.refractive_index.values,
    scatterer.diameter,
    values,
    shading="auto"
)
ax.set(
    xlabel="Permittivity",
    ylabel=r"Diameter [$\mu$m]",
    title="Coupling power [Watt]"
)
plt.colorbar(mappable=image)

plt.show()
