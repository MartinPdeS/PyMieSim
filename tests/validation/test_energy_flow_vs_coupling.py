import pytest
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single.scatterer import Sphere, CoreShell
from PyMieSim.single.source import Gaussian, PlaneWave
from PyMieSim.single.plottings import SystemPlotter

EPSILON0 = 8.854187817620389e-12 * ureg.farad / ureg.meter
C0 = 299792458.0 * ureg.meter / ureg.second

def test_energy_flow_vs_coupling_sphere():
    source = PlaneWave(
        wavelength=750 * ureg.nanometer,
        polarization=30 * ureg.degree,
        amplitude=1.0 * ureg.volt / ureg.meter,
    )

    scatterer = Sphere(
        diameter=550 * ureg.nanometer,
        source=source,
        refractive_index=1.5 * ureg.RIU,
        medium_refractive_index=1.0 * ureg.RIU,
    )

    # Define the detector (integrating sphere)
    detector = IntegratingSphere(
        sampling=5_000,
    )

    detector.print_properties(3)

    print(f"\n\n omega: {detector._cpp_mesh.omega}, sampling * d_omega: {detector._cpp_mesh.sampling * detector._cpp_mesh.d_omega}")

    # Calculate coupling and scattering efficiency (Qsca)
    coupling = detector.get_coupling(scatterer=scatterer)
    energy_flow = detector.get_energy_flow(scatterer)

    n = scatterer.medium_refractive_index

    intensity = (n * EPSILON0 * C0 * abs(source.amplitude)**2).to("watt/meter**2")

    scattered_power = (scatterer.Qsca * intensity * scatterer.cross_section).to("watt").to_compact()

    print("Qsca:", scatterer.Qsca)
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

def _test_energy_flow_vs_coupling_cylinder():
    source = Gaussian(
        wavelength=750 * ureg.nanometer,  # 750 nm
        polarization=30 * ureg.degree,  # Polarization in ureg.degrees
        optical_power=1 * ureg.watt,  # Power in ureg.watts
        NA=0.3 * ureg.AU,  # Numerical Aperture
    )

    # Define the scatterer (sphere)
    scatterer = CoreShell(
        core_diameter=1500 * ureg.nanometer,  # 1500 nm diameter
        shell_thickness=200 * ureg.nanometer,
        source=source,
        core_refractive_index=1.8 * ureg.RIU,  # Refractive index
        shell_refractive_index=1.5 * ureg.RIU,  # Refractive index
        medium_refractive_index=1.0 * ureg.RIU,  # Medium is water
    )

    # Define the detector (integrating sphere)
    detector = IntegratingSphere(
        sampling=5000,
    )

    # Calculate coupling and scattering efficiency (Qsca)
    coupling = detector.get_coupling(scatterer=scatterer)
    Qsca = scatterer.Qsca

    # Calculate energy flow
    energy_flow = detector.get_energy_flow(scatterer).to_compact()


    # Calculate scattered power
    scattered_power = (Qsca * source.peak_intensity * scatterer.cross_section).to_compact()

    # Check if the results are consistent
    assert np.isclose(
        coupling.to_compact(), scattered_power, atol=0, rtol=1e-2
    ), f"Mismatch betweend scattered power: {scattered_power} and coupling calculation: {coupling}"
    assert np.isclose(
        coupling.to_compact(), energy_flow, atol=0, rtol=1e-1
    ), f"Mismatch betweend energy flow: {energy_flow} and coupling calculation: {coupling}"

    # Assertions to validate the results
    assert coupling > 0, "Coupling should be a positive value."
    assert scattered_power > 0, "Scattered power should be a positive value."
    assert energy_flow > 0, "Energy flow should be a positive value."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
