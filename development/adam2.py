import numpy as np

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure
from PyOptik import UsualMaterial

source = Gaussian(
    wavelength=950e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

scatterer = Sphere(
    diameter=np.linspace(100e-9, 2000e-9, 20),
    index=1.4,
    medium_index=1,
    source=source
)

detector = CoherentMode(
    mode_number="HG11:00",
    NA=[0.2,0.3,0.4],
    phi_offset=np.linspace(0,360,200),
    gamma_offset=0,
    polarization_filter=None,
    sampling=300,
    rotation=0,  # Rotation of the mode field
)

experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

data = experiment.get(measure.coupling)

figure = data.plot(
    x=experiment.gamma_offset,
    std=experiment.polarization_filter
)

_ = figure.show()

"""
import numpy as np

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure
from PyOptik import UsualMaterial

source = Gaussian(
    wavelength=1310e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

scatterer = Sphere(
    diameter=500,
    index=1.4,
    medium_index=1,
    source=source
)

detector = CoherentMode(
    mode_number="HG11:00",
    NA=[0.2,0.3,0.4],
    phi_offset=np.linspace(0,360,200),
    gamma_offset=0,
    polarization_filter=None,
    sampling=500,
    rotation=0,  # Rotation of the mode field
)

experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

data = experiment.get(measure.coupling)

figure = data.plot(
    x=experiment.phi_offset,
    #std=experiment.diameter
)

_ = figure.show()"""