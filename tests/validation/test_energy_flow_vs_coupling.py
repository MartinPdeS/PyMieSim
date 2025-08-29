import pytest
import numpy as np
from TypedUnit import ureg
from PyOptik import Material

from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian


@pytest.fixture
def setup_simulation():
    # Define the source
    source = Gaussian(
        wavelength=750 * ureg.nanometer,  # 750 nm
        polarization=30 * ureg.degree,  # Polarization in ureg.degrees
        optical_power=1 * ureg.watt,  # Power in ureg.watts
        NA=0.3 * ureg.AU,  # Numerical Aperture
    )

    # Define the scatterer (sphere)
    scatterer = Sphere(
        diameter=1500 * ureg.nanometer,  # 1500 nm diameter
        source=source,
        property=1.4 * ureg.RIU,  # Refractive index
        medium_property=Material.water,  # Medium is water
    )

    # Define the detector (integrating sphere)
    detector = IntegratingSphere(
        sampling=1000,
    )

    return source, scatterer, detector


def test_simulation_results(setup_simulation):
    source, scatterer, detector = setup_simulation

    # Run simulation for Stokes parameters
    scatterer.get_stokes(distance=2 * ureg.meter, sampling=100)

    # Calculate coupling and scattering efficiency (Qsca)
    coupling = detector.get_coupling(scatterer=scatterer)
    Qsca = scatterer.Qsca

    # Calculate energy flow
    energy_flow = detector.get_energy_flow(scatterer, distance=1 * ureg.meter)

    # Calculate scattered power
    scattered_power = Qsca * source.peak_intensity * scatterer.cross_section

    # Check if the results are consistent
    assert np.isclose(
        coupling, scattered_power, atol=0, rtol=1e-2
    ), "Mismatch betweend scattered power: {scattered_power} and coupling calculation: {coupling}"
    assert np.isclose(
        coupling, energy_flow, atol=0, rtol=1e-1
    ), f"Mismatch betweend energy flow: {energy_flow} and coupling calculation: {coupling}"

    # Assertions to validate the results
    assert coupling > 0, "Coupling should be a positive value."
    assert scattered_power > 0, "Scattered power should be a positive value."
    assert energy_flow > 0, "Energy flow should be a positive value."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
