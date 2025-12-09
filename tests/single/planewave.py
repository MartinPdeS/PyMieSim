import pytest
from unittest.mock import patch
from TypedUnit import ureg

from PyMieSim.single.source import PlaneWave
from PyMieSim.single.polarization import Linear
from PyMieSim.single import plot_system


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


@patch("pyvista.Plotter.show")
def test_plotting(mock_show, source):
    """Test the plot method of PlaneWave to ensure it calls the show method once."""
    source.plot()
    mock_show.assert_called_once()


@patch("pyvista.Plotter.show")
def test_plot_system(mock_show):
    source = PlaneWave(
        amplitude=1 * ureg.volt / ureg.meter,
        polarization=0 * ureg.degree,
        wavelength=1550 * ureg.nanometer,
    )

    plot_system(source)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
