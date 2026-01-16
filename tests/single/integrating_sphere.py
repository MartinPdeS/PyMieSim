import pytest
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import IntegratingSphere

# Define sampling to be tested
samplings = [100, 200]


@pytest.fixture
def source():
    """Fixture to create a Gaussian source used across multiple tests."""
    return Gaussian(
        wavelength=750 * ureg.nanometer,  # Wavelength of the source in meters
        polarization=0 * ureg.degree,  # Polarization value
        optical_power=1 * ureg.watt,  # Optical power in watts
        NA=0.3 * ureg.AU,  # Numerical aperture
    )


@pytest.fixture
def setup_scatterer(source):
    """Fixture to create a scatterer with a provided source."""
    return Sphere(
        diameter=100 * ureg.nanometer,  # Diameter of the scatterer in meters
        source=source,  # Gaussian source from source fixture
        property=1.4 * ureg.RIU,  # Refractive index of the scatterer
        medium_property=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    )


@pytest.mark.parametrize("sampling", samplings)
def test_photodiode_with_sampling(sampling, setup_scatterer):
    """Test the Integrating Sphere detector with different sampling rates."""

    detector = IntegratingSphere(sampling=sampling)

    footprint = setup_scatterer.get_footprint(detector=detector)

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."

if __name__ == "__main__":
    pytest.main(["-W error", __file__])
