import pytest
from unittest.mock import patch
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import UnitPolarizationAngle



def test_gaussian_initialization():
    """Test the initialization of Gaussian source with different polarization inputs."""
    # Test with UnitPolarizationAngle
    gaussian1 = Gaussian(
        optical_power=1,
        NA=0.1,
        polarization=UnitPolarizationAngle(angle=0),
        wavelength=1550e-9
    )
    assert isinstance(gaussian1, Gaussian), "Failed to initialize Gaussian with UnitPolarizationAngle."

    # Test with scalar polarization
    gaussian2 = Gaussian(
        optical_power=1,
        NA=0.1,
        polarization=0,
        wavelength=1550e-9
    )
    assert isinstance(gaussian2, Gaussian), "Failed to initialize Gaussian with scalar polarization."




@pytest.fixture
def gaussian_source():
    """Fixture to create a Gaussian source with predefined parameters."""
    return Gaussian(
        optical_power=1,
        NA=0.1,
        polarization=0,
        wavelength=1550e-9
    )

@patch('pyvista.Plotter.show')
def test_gaussian_plotting(mock_show, gaussian_source):
    """Test the plot method of Gaussian to ensure it calls the show method once."""
    gaussian_source.plot()
    mock_show.assert_called_once()

if __name__ == "__main__":
    pytest.main([__file__])
