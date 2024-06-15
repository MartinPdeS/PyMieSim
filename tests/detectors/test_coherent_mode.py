import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode

# Define a list of mode numbers and rotation angles to be tested
mode_numbers = [
    "LP01", "LP11", "LP21",
    "LG01", "LG11", "LG21",
    "HG01", "HG11", "HG21"
]

rotations = [0, 90]


@pytest.fixture
def setup_source():
    """Fixture to create a Gaussian source used across multiple tests."""

    return Gaussian(
        wavelength=750e-9,  # Wavelength of the source in meters
        polarization=0,  # Polarization value
        optical_power=1,  # Optical power in watts
        NA=0.3  # Numerical aperture
    )


@pytest.fixture
def scatterer(setup_source):
    """Fixture to create a scatterer with a provided source."""

    return Sphere(
        diameter=100e-9,  # Diameter in meters
        source=setup_source,  # Source defined in the setup_source fixture
        index=1.4,  # Refractive index of the scatterer
        medium_index=1.0  # Refractive index of the medium
    )


@pytest.mark.parametrize('mode_number', mode_numbers)
@pytest.mark.parametrize('rotation', rotations)
def test_lp_modes(mode_number, rotation, scatterer):
    """Test different LP, LG, and HG modes with varying rotations."""

    detector = CoherentMode(
        mode_number=mode_number,
        NA=0.2,  # Numerical aperture for the detector
        sampling=100,  # Field sampling
        gamma_offset=0,  # Gamma offset
        phi_offset=0,  # Phi offset
        rotation=rotation  # Rotation angle
    )

    footprint = detector.get_footprint(scatterer=scatterer)

    assert footprint is not None, "Expected a valid footprint but got None."


if __name__ == "__main__":
    pytest.main()

# -
