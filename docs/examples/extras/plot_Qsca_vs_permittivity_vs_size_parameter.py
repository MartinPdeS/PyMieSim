"""
Scattering efficiency of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.
"""

import numpy
from PyMieSim.experiment import SphereSet, SourceSet, Setup
from PyMieSim import measure

from MPSPlots.render2D import SceneList


permitivity = numpy.linspace(-10, 50, 200)

index = numpy.sqrt(permitivity.astype(complex))

diameter = numpy.linspace(1e-9, 200e-9, 200)

source_set = SourceSet(
    wavelength=400e-9,
    polarization=90,
    amplitude=1
)


scatterer_set = SphereSet(
    diameter=diameter,
    index=index,
    n_medium=1
)

experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

data = experiment.Get(measures=measure.Qsca)

data = abs(data.array.squeeze())

figure = SceneList(unit_size=(6, 6))

ax = figure.append_ax(
    x_label="Permittivity",
    y_label=r'Diameter [$\mu$ m]',
    title="Scattering efficiency of a sphere"
)

_ = ax.add_mesh(
    x=permitivity,
    y=diameter,
    scalar=numpy.log(data),
    colormap='viridis',
    y_scale_factor=1e6,
    show_colorbar=True
)


_ = figure.show()
