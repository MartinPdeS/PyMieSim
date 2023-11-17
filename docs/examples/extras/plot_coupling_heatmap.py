"""
Coupling heatmap of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.
"""

import numpy
from PyMieSim.experiment import SphereSet, SourceSet, Setup, PhotodiodeSet
from PyMieSim import measure

from MPSPlots.render2D import SceneList


index = numpy.linspace(1.3, 2.1, 300)

diameter = numpy.linspace(1e-9, 2000e-9, 300)

source_set = SourceSet(
    wavelength=400e-9,
    linear_polarization=90,
    amplitude=1
)


scatterer_set = SphereSet(
    diameter=diameter,
    index=index,
    n_medium=1
)


detector_set = PhotodiodeSet(
    polarization_filter=None,
    NA=0.3,
    phi_offset=0,
    gamma_offset=0
)

experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=detector_set
)

data = experiment.Get(measure.coupling)

data = abs(data.array.squeeze())

figure = SceneList(unit_size=(6, 6))

ax = figure.append_ax(
    x_label="Permittivity",
    y_label=r'Diameter [$\mu$ m]',
    title="Coupling power [Watt]"
)

artist = ax.add_mesh(
    x=index,
    y=diameter,
    scalar=data,
    y_scale_factor=1e6,
)

_ = ax.add_colorbar(
    colormap='viridis',
    artist=artist,
    symmetric=False,
    numeric_format='%.3e'
)


_ = figure.show()
