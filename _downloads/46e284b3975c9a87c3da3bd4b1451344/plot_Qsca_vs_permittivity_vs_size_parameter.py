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

from PyMieSim import measure

from MPSPlots.render2D import SceneList


permitivity = numpy.linspace(-10, 50, 400)

index = numpy.sqrt(permitivity.astype(complex))

diameter = numpy.linspace(1e-9, 200e-9, 400)

source_set = Gaussian(
    wavelength=400e-9,
    polarization_value=90,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)


scatterer_set = Sphere(
    diameter=diameter,
    index=index,
    n_medium=1,
    source_set=source_set
)

experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

data = experiment.get(measure.Qsca)

data = abs(data.y.values.squeeze())

figure = SceneList(unit_size=(6, 6))

ax = figure.append_ax(
    x_label="Permittivity",
    y_label=r'Diameter [$\mu$ m]',
    title="Scattering efficiency of a sphere"
)

artist = ax.add_mesh(
    x=permitivity,
    y=diameter,
    scalar=numpy.log(data),
    y_scale_factor=1e6,
)

_ = ax.add_colorbar(artist=artist, colormap='viridis')

_ = figure.show()
