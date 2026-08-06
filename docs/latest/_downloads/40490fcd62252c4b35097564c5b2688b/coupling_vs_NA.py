"""
Sphere: Coupling vs numerical aperture
======================================
"""

import numpy as np
from PyMieSim import (
    Experiment,
    Gaussian,
    GaussianSet,
    PhotodiodeSet,
    PolarizationSet,
    PolarizationState,
    RightCircular,
    SellmeierMaterial,
    Simulation,
    Sphere,
    SphereSet,
    ureg,
)

source = GaussianSet(
    wavelength=[500] * ureg.nanometer,
    polarization=PolarizationSet(angles=[0] * ureg.degree),
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = SphereSet(
    diameter=[500e-9] * ureg.meter,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = PhotodiodeSet(
    numerical_aperture=np.linspace(0.1, 1, 150),
    phi_offset=[0] * ureg.degree,
    gamma_offset=[0, 10] * ureg.degree,
    sampling=[2000]
)

setup = Experiment(
    scatterer_set=scatterer,
    source_set=source,
    detector_set=detector
)

result = setup.get("coupling", drop_unique_level=True)

result.plot(x="detector:NA")

single_source = Gaussian(
    wavelength=950 * ureg.nanometer,
    polarization=PolarizationState(angle=0 * ureg.degree),
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2,
)

single_scatterer = Sphere(
    diameter=500 * ureg.nanometer,
    material=1.5,
    medium=1,
)

setup = Simulation(
    source=single_source,
    scatterer=single_scatterer
)

print(setup.get("Qsca"))
