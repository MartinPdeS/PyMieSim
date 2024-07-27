"""
Coupling heatmap of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.
"""

import numpy

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure

from MPSPlots.render2D import SceneList


index = numpy.linspace(1.3, 2.1, 300)

diameter = numpy.linspace(1e-9, 2000e-9, 200)

source = Gaussian(
    wavelength=400e-9,
    polarization=90,
    optical_power=1e-3,
    NA=0.2
)


scatterer = Sphere(
    diameter=diameter,
    index=index,
    medium_index=1,
    source=source
)


detector = Photodiode(
    polarization_filter=None,
    NA=0.3,
    phi_offset=0,
    gamma_offset=0,
    sampling=400
)

experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

data = experiment.get(measure.coupling)

data = abs(data.y.values.squeeze())

figure = SceneList(unit_size=(6, 6))

ax = figure.append_ax(
    x_label="Permittivity",
    y_label=r'Diameter [$\mu$m]',
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
