"""
InfiniteCylinder: Goniometer
============================

This example demonstrates how to use a goniometer setup to measure and visualize the coupling efficiency as a function of angular displacement for cylindrical scatterers using PyMieSim.
"""

import PyOptik
from PyMieSim import from_pyoptik
import numpy as np
from PyMieSim import (
    Measure,
    ureg,
    PhotodiodeSet,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)



polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=[2000] * ureg.nanometer,
    material=[from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7"))],
    medium=[1],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.5, 0.3, 0.1, 0.05],
    phi_offset=np.linspace(-180, 180, 200) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[400],
    polarization_filter=None,
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.get(Measure.COUPLING)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.attrs["units"])

result.plot(x="detector:phi_offset")
