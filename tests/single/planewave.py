import pytest
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.source import PlaneWave
from PyMieSim.single.polarization import Linear


def test_planewave_initialization():
    """Test the initialization of PlaneWave source with different polarization inputs."""
    # Test with Linear
    planewave_1 = PlaneWave(
        amplitude=1 * ureg.volt / ureg.meter,
        polarization=Linear(element=0 * ureg.degree),
        wavelength=1550 * ureg.nanometer,
    )
    assert isinstance(
        planewave_1, PlaneWave
    ), "Failed to initialize PlaneWave with Linear."

    # Test with scalar polarization
    planewave_2 = PlaneWave(
        amplitude=1 * ureg.volt / ureg.meter,
        polarization=0 * ureg.degree,
        wavelength=1550 * ureg.nanometer,
    )
    assert isinstance(
        planewave_2, PlaneWave
    ), "Failed to initialize PlaneWave with scalar polarization."


@pytest.fixture
def source():
    """Fixture to create a PlaneWave source with predefined parameters."""
    return PlaneWave(
        amplitude=1 * ureg.volt / ureg.meter,
        polarization=0 * ureg.degree,
        wavelength=1550 * ureg.nanometer,
    )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
