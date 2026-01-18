import pytest
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single.scatterer import Sphere, CoreShell
from PyMieSim.single.source import Gaussian


def test_energy_flow_vs_coupling_sphere():
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
        property=1.8 * ureg.RIU,  # Refractive index
        medium_property=1.4 * ureg.RIU,  # Medium is water
    )

    # Define the detector (integrating sphere)
    detector = IntegratingSphere(
        sampling=5000,
    )

    # Run simulation for Stokes parameters
    scatterer.get_stokes(distance=2 * ureg.meter, sampling=100)

    # Calculate coupling and scattering efficiency (Qsca)
    coupling = detector.get_coupling(scatterer=scatterer)
    Qsca = scatterer.Qsca

    # Calculate energy flow
    energy_flow = detector.get_energy_flow(scatterer)


    # Calculate scattered power
    scattered_power = Qsca * source.peak_intensity * scatterer.cross_section

    # Check if the results are consistent
    assert np.isclose(
        coupling.to_compact(), scattered_power.to_compact(), atol=0, rtol=1e-2
    ), "Mismatch betweend scattered power: {scattered_power} and coupling calculation: {coupling}"
    assert np.isclose(
        coupling.to_compact(), energy_flow.to_compact(), atol=0, rtol=1e-1
    ), f"Mismatch betweend energy flow: {energy_flow} and coupling calculation: {coupling}"

    # Assertions to validate the results
    assert coupling > 0, "Coupling should be a positive value."
    assert scattered_power > 0, "Scattered power should be a positive value."
    assert energy_flow > 0, "Energy flow should be a positive value."

def test_energy_flow_vs_coupling_cylinder():
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
        core_property=1.8 * ureg.RIU,  # Refractive index
        shell_property=1.5 * ureg.RIU,  # Refractive index
        medium_property=1.4 * ureg.RIU,  # Medium is water
    )

    # Define the detector (integrating sphere)
    detector = IntegratingSphere(
        sampling=5000,
    )

    # Run simulation for Stokes parameters
    scatterer.get_stokes(distance=2 * ureg.meter, sampling=100)

    # Calculate coupling and scattering efficiency (Qsca)
    coupling = detector.get_coupling(scatterer=scatterer)
    Qsca = scatterer.Qsca

    # Calculate energy flow
    energy_flow = detector.get_energy_flow(scatterer)


    # Calculate scattered power
    scattered_power = Qsca * source.peak_intensity * scatterer.cross_section

    # Check if the results are consistent
    assert np.isclose(
        coupling.to_compact(), scattered_power.to_compact(), atol=0, rtol=1e-2
    ), f"Mismatch betweend scattered power: {scattered_power} and coupling calculation: {coupling}"
    assert np.isclose(
        coupling.to_compact(), energy_flow.to_compact(), atol=0, rtol=1e-1
    ), f"Mismatch betweend energy flow: {energy_flow} and coupling calculation: {coupling}"

    # Assertions to validate the results
    assert coupling > 0, "Coupling should be a positive value."
    assert scattered_power > 0, "Scattered power should be a positive value."
    assert energy_flow > 0, "Energy flow should be a positive value."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
