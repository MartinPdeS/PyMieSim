from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment import Setup, measure

from PyOptik import UsualMaterial
import numpy

source = Gaussian(
    wavelength=[400e-9, 488e-9, 638e-9],
    polarization_value=numpy.linspace(0, 180, 15),
    polarization_type='linear',
    optical_power=1,
    NA=0.3
)

sphere = Sphere(
    diameter=numpy.linspace(50e-9, 400e-9),
    material=UsualMaterial.Polystyrene,
    source=source,
    medium_material=UsualMaterial.Water
)

detector = Photodiode(
    sampling=500,
    NA=0.2,
    gamma_offset=0,
    phi_offset=[0, 45, 90],
    polarization_filter=0,
)


experiment = Setup(
    source=source,
    scatterer=sphere,
    detector=detector
)

data = experiment.get(measure.coupling)


figure = data.plot(
    x=experiment.diameter,
    std=experiment.polarization_value
)

figure.show()
