import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode

# Define different sampling for testing
samplings = [100, 200, 300, 400]


@pytest.fixture
def setup_source():
    """Fixture to create a Gaussian source that can be reused in different tests."""

    return Gaussian(
        wavelength=750e-9,  # Wavelength of the source in meters
        polarization=0,  # Polarization value
        optical_power=1,  # Optical power in watts
        NA=0.3  # Numerical aperture
    )


@pytest.fixture
def setup_scatterer(setup_source):
    """Fixture to create a scatterer using the Gaussian source defined above."""

    return Sphere(
        diameter=100e-9,  # Diameter of the scatterer in meters
        source=setup_source,  # Source from setup_source fixture
        index=1.4,  # Refractive index of the scatterer
        medium_index=1.0  # Refractive index of the surrounding medium
    )


@pytest.mark.parametrize('sampling', samplings)
def test_photodiode_sampling(sampling, setup_scatterer):
    """Test the Photodiode detector with various sampling rates."""

    detector = Photodiode(
        NA=0.2,  # Numerical aperture of the detector
        sampling=sampling,  # Field sampling
        gamma_offset=0,  # Gamma offset
        phi_offset=0  # Phi offset
    )

    # Perform the operation to be tested
    footprint = detector.get_footprint(scatterer=setup_scatterer)

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."


if __name__ == "__main__":
    pytest.main()


# -
