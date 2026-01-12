import pytest
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import Linear
from PyMieSim.single import plot_system


def test_gaussian_initialization():
    """Test the initialization of Gaussian source with different polarization inputs."""
    # Test with UnitPolarizationAngle
    gaussian1 = Gaussian(
        optical_power=1 * ureg.watt,
        NA=0.1 * ureg.AU,
        polarization=Linear(element=0 * ureg.degree),
        wavelength=1550 * ureg.nanometer,
    )
    assert isinstance(
        gaussian1, Gaussian
    ), "Failed to initialize Gaussian with UnitPolarizationAngle."

    # Test with scalar polarization
    gaussian2 = Gaussian(
        optical_power=1 * ureg.watt,
        NA=0.1 * ureg.AU,
        polarization=0 * ureg.degree,
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
        NA=0.1 * ureg.AU,
        polarization=0 * ureg.degree,
        wavelength=1550 * ureg.nanometer,
    )


@patch("pyvista.Plotter.show")
def test_gaussian_plotting(mock_show, gaussian_source):
    """Test the plot method of Gaussian to ensure it calls the show method once."""
    gaussian_source.plot()
    mock_show.assert_called_once()


@patch("pyvista.Plotter.show")
def test_plot_system(mock_show):
    source = Gaussian(
        optical_power=1 * ureg.watt,
        NA=0.1 * ureg.AU,
        polarization=0 * ureg.degree,
        wavelength=1550 * ureg.nanometer,
    )

    plot_system(source)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
