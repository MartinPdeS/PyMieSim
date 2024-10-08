import pytest
import numpy as np
from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.units import AU, degree, nanometer, watt, RIU, meter
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyOptik import Material

@pytest.fixture
def setup_simulation():
    # Define the source
    source = Gaussian(
        wavelength=750 * nanometer,  # 750 nm
        polarization=30 * degree,  # Polarization in degrees
        optical_power=1 * watt,  # Power in watts
        NA=0.3 * AU  # Numerical Aperture
    )

    # Define the scatterer (sphere)
    scatterer = Sphere(
        diameter=1500 * nanometer,  # 1500 nm diameter
        source=source,
        property=1.4 * RIU,  # Refractive index
        medium_property=Material.water  # Medium is water
    )

    # Define the detector (integrating sphere)
    detector = IntegratingSphere(
        sampling=1000 * AU,
        polarization_filter=None
    )

    return source, scatterer, detector

def test_simulation_results(setup_simulation):
    source, scatterer, detector = setup_simulation

    # Run simulation for Stokes parameters
    scatterer.get_stokes(distance=2 * meter, sampling=100)

    # Calculate coupling and scattering efficiency (Qsca)
    coupling = detector.coupling(scatterer=scatterer)
    Qsca = scatterer.Qsca

    # Calculate energy flow
    energy_flow = detector.get_energy_flow(scatterer, distance=1 * meter)

    # Calculate scattered power
    scattered_power = Qsca * source.peak_intensity * scatterer.cross_section

    # Print the results
    print(coupling, scattered_power, energy_flow)

    # Assertions to validate the results
    assert coupling > 0, "Coupling should be a positive value."
    assert scattered_power > 0, "Scattered power should be a positive value."
    assert energy_flow > 0, "Energy flow should be a positive value."

if __name__ == '__main__':
    pytest.main([__file__])