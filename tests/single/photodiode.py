import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import Setup


@pytest.fixture
def source():
    """Fixture to create a Gaussian source that can be reused in different tests."""
    return Gaussian(
        wavelength=750 * ureg.nanometer,  # Wavelength of the source in meters
        polarization=PolarizationState(angle=0 * ureg.degree),  # Polarization value
        optical_power=1 * ureg.watt,  # Optical power in ureg.watts
        numerical_aperture=0.3 * ureg.AU,  # Numerical aperture
    )


@pytest.fixture
def setup_scatterer():
    """Fixture to create a scatterer using the Gaussian source defined above."""
    return Sphere(
        diameter=100 * ureg.nanometer,  # Diameter of the scatterer in meters
        material=1.4 * ureg.RIU,  # Refractive index of the scatterer
        medium=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    )


@pytest.fixture
def photodiode():
    """Test the Photodiode detector with various sampling rates."""
    return Photodiode(
        numerical_aperture=0.2 * ureg.AU,  # Numerical aperture of the detector
        sampling=30,  # Field sampling
        gamma_offset=0 * ureg.degree,  # Gamma offset
        phi_offset=0 * ureg.degree,  # Phi offset
        medium=1.0 * ureg.RIU
    )


def test_photodiode_sampling(photodiode, source, setup_scatterer):
    """Test the Photodiode detector with various sampling rates."""

    setup = Setup(
        scatterer=setup_scatterer,
        source=source,
        detector=photodiode,
    )
    # Perform the operation to be tested
    footprint = setup.get_representation("footprint")

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."


def test_fails_initialization():
    with pytest.raises(Exception):
        Photodiode(
            numerical_aperture=0.2 * ureg.AU,
            block_numerical_aperture=0.3 * ureg.AU,
            gamma_offset=0 * ureg.degree,
            phi_offset=0 * ureg.degree,
            medium=1.0 * ureg.RIU
        )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
