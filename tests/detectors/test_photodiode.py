import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from unittest.mock import patch
from PyMieSim.units import nanometer, degree, watt, AU, RIU

@pytest.fixture
def setup_source():
    """Fixture to create a Gaussian source that can be reused in different tests."""

    return Gaussian(
        wavelength=750 * nanometer,  # Wavelength of the source in meters
        polarization=0 * degree,  # Polarization value
        optical_power=1 * watt,  # Optical power in watts
        NA=0.3 * AU # Numerical aperture
    )


@pytest.fixture
def setup_scatterer(setup_source):
    """Fixture to create a scatterer using the Gaussian source defined above."""

    return Sphere(
        diameter=100 * nanometer,  # Diameter of the scatterer in meters
        source=setup_source,  # Source from setup_source fixture
        property=1.4 * RIU,  # Refractive index of the scatterer
        medium_property=1.0 * RIU # Refractive index of the surrounding medium
    )

@pytest.fixture
def photodiode():
    """Test the Photodiode detector with various sampling rates."""

    return Photodiode(
        NA=0.2 * AU,  # Numerical aperture of the detector
        sampling=30 * AU,  # Field sampling
        gamma_offset=0 * degree,  # Gamma offset
        phi_offset=0 * degree # Phi offset
    )


@patch('pyvista.Plotter.show')
def test_photodiode_plot(mock_show, photodiode):
    photodiode.plot()


def test_photodiode_sampling(photodiode, setup_scatterer):
    """Test the Photodiode detector with various sampling rates."""

    # Perform the operation to be tested
    footprint = photodiode.get_footprint(scatterer=setup_scatterer)

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."


if __name__ == "__main__":
    pytest.main([__file__])


# -
