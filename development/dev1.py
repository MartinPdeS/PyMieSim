import pytest
import numpy as np
from PyMieSim.units import ureg
from PyOptik import Material

from PyMieSim.single.detector import IntegratingSphere, Photodiode
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian



source = Gaussian(
    wavelength=750 * ureg.nanometer,  # 750 nm
    polarization=30 * ureg.degree,  # Polarization in ureg.degrees
    optical_power=1 * ureg.watt,  # Power in ureg.watts
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

scatterer = Sphere(
    diameter=1500 * ureg.nanometer,
    source=source,
    property=1.8 * ureg.RIU,
    medium_property=1.01 * ureg.RIU
)

# Define the detector (integrating sphere)
detector = IntegratingSphere(
    sampling=5000,
)

from PyMieSim.single import SystemPlotter
plotter = SystemPlotter()


# plotter.plot(detector)

coupling = detector.get_coupling(scatterer=scatterer)
print("coupling", coupling)
Qsca = scatterer.Qsca

energy_flow = detector.get_energy_flow(scatterer)


# Calculate scattered power
scattered_power = Qsca * source.peak_intensity * scatterer.cross_section
print('0-------------', scattered_power.to_compact() / coupling.to_compact())

# Check if the results are consistent
assert np.isclose(
    coupling.to_compact(), scattered_power.to_compact(), atol=0, rtol=1e-2
), f"Mismatch betweend scattered power: {scattered_power.to_compact()} and coupling calculation: {coupling.to_compact()}"
assert np.isclose(
    coupling.to_compact(), energy_flow.to_compact(), atol=0, rtol=1e-1
), f"Mismatch betweend energy flow: {energy_flow.to_compact()} and coupling calculation: {coupling.to_compact()}"

# Assertions to validate the results
assert coupling > 0, "Coupling should be a positive value."
assert scattered_power > 0, "Scattered power should be a positive value."
assert energy_flow > 0, "Energy flow should be a positive value."
