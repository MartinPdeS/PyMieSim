import numpy as np
from PyMieSim.experiment import Sphere, Gaussian, Photodiode, Setup
from PyMieSim import measure

source = Gaussian(
    wavelength=400e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.3
)

scatterer_set = Sphere(
    diameter=800e-9,
    index=1.44,
    medium_index=1,
    source=source
)

detector_set = Photodiode(
    NA=[0.1, 0.2],
    phi_offset=np.linspace(-180, 180, 200),
    gamma_offset=0,
    sampling=[100, 300],
    polarization_filter=None
)

setup = Setup(
    source=source,
    scatterer=scatterer_set,
    detector=detector_set
)

data_set = setup.get(measure.coupling)

figure = data_set.plot(
    x=setup.phi_offset,
)


_ = figure.show()