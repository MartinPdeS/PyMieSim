"""
Inspect scalar scattering properties
====================================

Print scalar efficiencies, cross sections, and the asymmetry factor for one
isolated sphere. Efficiencies and asymmetry are dimensionless; cross sections
carry area units. The real refractive index gives negligible absorption, so
``Qext`` and ``Qsca`` agree.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim import (
    Measure,
    ureg,
    Gaussian,
    PolarizationState,
    Sphere,
    Simulation,
)

polarization = PolarizationState(angle=30 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = Sphere(
    diameter=1500 * ureg.nanometer,
    material=1.4,
    medium=1.0,
)

setup = Simulation(
    scatterer=scatterer,
    source=source,
)

for measure in (
    Measure.QEXT, Measure.QSCA, Measure.QABS, Measure.QBACK, Measure.G,
    Measure.CEXT, Measure.CSCA, Measure.CABS, Measure.CBACK,
):
    value = setup.run(measure)
    print(f"{measure.value}: {value}")


