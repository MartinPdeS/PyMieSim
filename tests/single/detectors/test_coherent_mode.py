import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode
from unittest.mock import patch
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyMieSim.single import plot_system


# Define a list of mode numbers and rotation angles to be tested
mode_numbers = [
    "LP01", "LP11", "LP21",
    "LG01", "LG11", "LG21",
    "HG01", "HG11", "HG21"
]


@pytest.fixture
def source():
    """Fixture to create a Gaussian source used across multiple tests."""

    return Gaussian(
        wavelength=750 * nanometer,  # Wavelength of the source in meters
        polarization=0 * degree,  # Polarization value
        optical_power=1 * watt,  # Optical power in watts
        NA=0.3 * AU  # Numerical aperture
    )


@pytest.fixture
def scatterer(source):
    """Fixture to create a scatterer with a provided source."""

    return Sphere(
        diameter=100 * nanometer,  # Diameter in meters
        source=source,  # Source defined in the setup_source fixture
        property=1.4 * RIU,  # Refractive index of the scatterer
        medium_property=1.0 * RIU  # Refractive index of the medium
    )


@pytest.mark.parametrize('mode_number', mode_numbers)
def test_lp_modes(mode_number, scatterer):
    """Test different LP, LG, and HG modes with varying rotations."""

    detector = CoherentMode(
        mode_number=mode_number,
        NA=0.2 * AU,  # Numerical aperture for the detector
        sampling=100,  # Field sampling
        gamma_offset=0 * degree,  # Gamma offset
        phi_offset=0 * degree,  # Phi offset
        rotation=0 * degree  # Rotation angle
    )

    footprint = detector.get_footprint(scatterer=scatterer)

    assert footprint is not None, "Expected a valid footprint but got None."


@patch('pyvista.Plotter.show')
def test_plot_system(mock_show):
    detector = CoherentMode(
        mode_number="LP01",
        NA=0.2 * AU,  # Numerical aperture for the detector
        sampling=100,  # Field sampling
        gamma_offset=0 * degree,  # Gamma offset
        phi_offset=0 * degree,  # Phi offset
        rotation=0 * degree  # Rotation angle
    )


    plot_system(detector)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
