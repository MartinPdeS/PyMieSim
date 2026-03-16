"""
Sphere: Coupling vs numerical aperture
======================================
"""

import numpy as np
from PyMieSim.units import ureg
from PyMieSim.material import SellmeierMaterial

from PyMieSim import experiment
from PyMieSim import single
from PyMieSim.polarization import PolarizationState

source = experiment.source_set.GaussianSet(
    wavelength=[500] * ureg.nanometer,
    polarization=experiment.polarization_set.PolarizationSet(angles=[0] * ureg.degree),
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = experiment.scatterer_set.SphereSet(
    diameter=[500e-9] * ureg.meter,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = experiment.detector_set.PhotodiodeSet(
    numerical_aperture=np.linspace(0.1, 1, 150),
    phi_offset=[0] * ureg.degree,
    gamma_offset=[0, 10] * ureg.degree,
    sampling=[2000]
)

setup = experiment.Setup(
    scatterer_set=scatterer,
    source_set=source,
    detector_set=detector
)

dataframe = setup.get("coupling", drop_unique_level=True)

dataframe.plot(x="detector:NA")

single_source = single.source.Gaussian(
    wavelength=950 * ureg.nanometer,
    polarization=PolarizationState(angle=0 * ureg.degree),
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2,
)

single_scatterer = single.scatterer.Sphere(
    diameter=500 * ureg.nanometer,
    material=1.5,
    medium=1,
)

setup = single.setup.Setup(
    source=single_source,
    scatterer=single_scatterer
)

print(setup.get("Qsca"))
