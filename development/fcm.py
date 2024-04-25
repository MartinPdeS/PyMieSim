import numpy
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import UsualMaterial
from PyMieSim import measure

source_set = Gaussian(
    wavelength=600e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

scatterer_set = Sphere(
    diameter=3e-6,
    index=1.5,
    medium_material=UsualMaterial.Water,
    source_set=source_set
)

detector_set = CoherentMode(
    mode_number='LP11',
    NA=[0.1],
    phi_offset=0,
    gamma_offset=0,
    sampling=800,
    rotation=numpy.linspace(0, 180, 100),
    polarization_filter=None
)

experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=detector_set
)

data = experiment.get(measure.coupling)

figure = data.plot(
    x=experiment.rotation,
    # y_scale='log',
    # std=experiment.index
)

_ = figure.show()
