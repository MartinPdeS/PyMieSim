import pytest
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import PlaneWave
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.setup import Setup

EPSILON0 = 8.854187817620389e-12 * ureg.farad / ureg.meter
C0 = 299792458.0 * ureg.meter / ureg.second

def test_energy_flow_vs_coupling_sphere():
    source = PlaneWave(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=30 * ureg.degree),
        amplitude=1.0 * ureg.volt / ureg.meter,
    )

    scatterer = Sphere(
        diameter=550 * ureg.nanometer,
        material=1.5,
        medium=1.0,
    )

    detector = IntegratingSphere(
        sampling=5_000,
    )

    setup = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector
    )

    # Calculate coupling and scattering efficiency (Qsca)
    coupling = setup.get("coupling")

    energy_flow = detector.get_energy_flow(scatterer=scatterer, source=source)

    n = scatterer.medium.refractive_index

    intensity = (0.5 * n * EPSILON0 * C0 * abs(source.amplitude)**2).to("watt/meter**2")

    scattered_power = (setup.get("Qsca") * intensity * setup.get("cross_section")).to("watt").to_compact()

    print("Qsca:", setup.get("Qsca"))
    print("coupling:", coupling.to_compact())
    print("energy_flow:", energy_flow.to_compact())
    print("theory:", scattered_power)
    print("theory/energy:", (scattered_power / energy_flow).to("dimensionless"))

    print(
        "Coupling: ", coupling.to_compact(),
        "\nScattered power: ", scattered_power,
        "\nEnergy flow: ", energy_flow.to_compact(),
        "\nscattered_power / Energy flow: ", (scattered_power / energy_flow).to_compact(),
        "\nEnergy flow / scattered_power: ", (energy_flow / scattered_power).to_compact()
    )
    # Check if the results are consistent
    assert np.isclose(
        coupling, scattered_power, atol=0, rtol=1e-2
    ), f"Mismatch betweend scattered power: {scattered_power} and coupling calculation: {coupling}"
    assert np.isclose(
        coupling, energy_flow, atol=0, rtol=1e-1
    ), f"Mismatch betweend energy flow: {energy_flow} and coupling calculation: {coupling}"

    # Assertions to validate the results
    assert coupling > 0, "Coupling should be a positive value."
    assert scattered_power > 0, "Scattered power should be a positive value."
    assert energy_flow > 0, "Energy flow should be a positive value."

if __name__ == "__main__":
    pytest.main(["-W error", __file__])
