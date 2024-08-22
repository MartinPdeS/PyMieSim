import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import IntegratingSphere

# Define sampling to be tested
samplings = [100, 200]


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
def setup_scatterer(setup_source):
    """Fixture to create a scatterer with a provided source."""

    return Sphere(
        diameter=100e-9,  # Diameter of the scatterer in meters
        source=setup_source,  # Gaussian source from setup_source fixture
        index=1.4,  # Refractive index of the scatterer
        medium_index=1.0  # Refractive index of the surrounding medium
    )


@pytest.mark.parametrize('sampling', samplings)
def test_photodiode_with_sampling(sampling, setup_scatterer):
    """Test the Integrating Sphere detector with different sampling rates."""

    detector = IntegratingSphere(sampling=sampling)

    footprint = detector.get_footprint(scatterer=setup_scatterer)

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."


if __name__ == "__main__":
    pytest.main()

# -
