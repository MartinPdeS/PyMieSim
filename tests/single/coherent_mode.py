import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import Setup

# Define a list of mode numbers and rotation angles to be tested
mode_numbers = ["LP01", "LP11", "LP21", "LG01", "LG11", "LG21", "HG01", "HG11", "HG21"]


@pytest.mark.parametrize("mode_number", mode_numbers)
def test_lp_modes(mode_number):
    """Test different LP, LG, and HG modes with varying rotations."""

    source = Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        material=1.4 * ureg.RIU,
        medium=1.0 * ureg.RIU,
    )

    detector = CoherentMode(
        mode_number=mode_number,
        numerical_aperture=0.2 * ureg.AU,
        sampling=100,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        rotation=0 * ureg.degree,
    )

    setup = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector,
    )

    footprint = setup.get_representation("footprint")

    assert footprint is not None, "Expected a valid footprint but got None."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
