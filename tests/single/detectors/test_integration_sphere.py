import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# Define sampling to be tested
samplings = [100, 200] * AU


@pytest.fixture
def setup_source():
    """Fixture to create a Gaussian source used across multiple tests."""

    return Gaussian(
        wavelength=750 * nanometer,  # Wavelength of the source in meters
        polarization=0 * degree,  # Polarization value
        optical_power=1 * watt,  # Optical power in watts
        NA=0.3 * AU  # Numerical aperture
    )


@pytest.fixture
def setup_scatterer(setup_source):
    """Fixture to create a scatterer with a provided source."""

    return Sphere(
        diameter=100 * nanometer,  # Diameter of the scatterer in meters
        source=setup_source,  # Gaussian source from setup_source fixture
        property=1.4 * RIU,  # Refractive index of the scatterer
        medium_property=1.0 * RIU  # Refractive index of the surrounding medium
    )


@pytest.mark.parametrize('sampling', samplings)
def test_photodiode_with_sampling(sampling, setup_scatterer):
    """Test the Integrating Sphere detector with different sampling rates."""

    detector = IntegratingSphere(sampling=sampling)

    footprint = detector.get_footprint(scatterer=setup_scatterer)

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
