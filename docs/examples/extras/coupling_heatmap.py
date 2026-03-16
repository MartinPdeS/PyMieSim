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

from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import SellmeierMedium

polarization_state = PolarizationSet(
    angles=[90] * ureg.degree,  # Linear polarization at 90 degrees
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_state,  # Linear polarization at 90 degrees
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

material_values = numpy.linspace(1.3, 2.1, 100)

scatterer = SphereSet(
    diameter=numpy.linspace(1, 2000, 100) * ureg.nanometer,
    material=material_values,
    medium=[SellmeierMedium("water")],
)

detector = PhotodiodeSet(
    polarization_filter=[0] * ureg.degree,
    numerical_aperture=[0.3],
    phi_offset=[0] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[400]
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

values = experiment.get("coupling", add_units=False, as_numpy=True)


figure, ax = plt.subplots(1, 1)

image = ax.pcolormesh(
    material_values,
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
