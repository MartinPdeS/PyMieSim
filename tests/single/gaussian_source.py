import pytest
from PyMieSim.units import ureg

from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState


def test_gaussian_initialization():
    """Test the initialization of Gaussian source with different polarization inputs."""
    # Test with UnitPolarizationAngle
    gaussian1 = Gaussian(
        optical_power=1 * ureg.watt,
        numerical_aperture=0.1 * ureg.AU,
        polarization=PolarizationState(angle=0 * ureg.degree),
        wavelength=1550 * ureg.nanometer,
    )
    assert isinstance(
        gaussian1, Gaussian
    ), "Failed to initialize Gaussian with UnitPolarizationAngle."

    # Test with scalar polarization
    gaussian2 = Gaussian(
        optical_power=1 * ureg.watt,
        numerical_aperture=0.1 * ureg.AU,
        polarization=PolarizationState(angle=0 * ureg.degree),
        wavelength=1550 * ureg.nanometer,
    )
    assert isinstance(
        gaussian2, Gaussian
    ), "Failed to initialize Gaussian with scalar polarization."


@pytest.fixture
def gaussian_source():
    """Fixture to create a Gaussian source with predefined parameters."""
    return Gaussian(
        optical_power=1 * ureg.watt,
        numerical_aperture=0.1 * ureg.AU,
        polarization=PolarizationState(angle=0 * ureg.degree),
        wavelength=1550 * ureg.nanometer,
    )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
