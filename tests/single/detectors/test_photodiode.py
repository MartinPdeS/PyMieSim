import pytest
from TypedUnit import ureg
from unittest.mock import patch

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import plot_system


@pytest.fixture
def source():
    """Fixture to create a Gaussian source that can be reused in different tests."""
    return Gaussian(
        wavelength=750 * ureg.nanometer,  # Wavelength of the source in meters
        polarization=0 * ureg.degree,  # Polarization value
        optical_power=1 * ureg.watt,  # Optical power in ureg.watts
        NA=0.3 * ureg.AU,  # Numerical aperture
    )


@pytest.fixture
def setup_scatterer(source):
    """Fixture to create a scatterer using the Gaussian source defined above."""
    return Sphere(
        diameter=100 * ureg.nanometer,  # Diameter of the scatterer in meters
        source=source,  # Source from source fixture
        property=1.4 * ureg.RIU,  # Refractive index of the scatterer
        medium_property=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    )


@pytest.fixture
def photodiode():
    """Test the Photodiode detector with various sampling rates."""
    return Photodiode(
        NA=0.2 * ureg.AU,  # Numerical aperture of the detector
        sampling=30,  # Field sampling
        gamma_offset=0 * ureg.degree,  # Gamma offset
        phi_offset=0 * ureg.degree,  # Phi offset
    )


def test_photodiode_sampling(photodiode, setup_scatterer):
    """Test the Photodiode detector with various sampling rates."""

    # Perform the operation to be tested
    footprint = photodiode.get_footprint(scatterer=setup_scatterer)

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."


def test_fails_initialization():
    with pytest.raises(Exception):
        Photodiode(
            NA=0.2 * ureg.AU,
            block_NA=0.3 * ureg.AU,
            gamma_offset=0 * ureg.degree,
            phi_offset=0 * ureg.degree,
        )


@patch("pyvista.Plotter.show")
def test_plot_system(mock_show):
    detector = Photodiode(
        NA=0.2 * ureg.AU,  # Numerical aperture of the detector
        sampling=30,  # Field sampling
        gamma_offset=0 * ureg.degree,  # Gamma offset
        phi_offset=0 * ureg.degree,  # Phi offset
    )

    plot_system(detector)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
