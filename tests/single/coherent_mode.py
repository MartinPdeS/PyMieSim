import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single.representations import Footprint

# Define a list of mode numbers and rotation angles to be tested
mode_numbers = ["LP01", "LP11", "LP21", "LG01", "LG11", "LG21", "HG01", "HG11", "HG21"]


@pytest.mark.parametrize("mode_number", mode_numbers)
def test_lp_modes(mode_number):
    """Test different LP, LG, and HG modes with varying rotations."""

    source = Gaussian(
        wavelength=750 * ureg.nanometer,  # Wavelength of the source in meters
        polarization=0 * ureg.degree,  # Polarization value
        optical_power=1 * ureg.watt,  # Optical power in watts
        NA=0.3 * ureg.AU,  # Numerical aperture
    )

    scatterer = Sphere(
        diameter=100 * ureg.nanometer,  # Diameter in meters
        source=source,  # Source defined in the setup_source fixture
        refractive_index=1.4 * ureg.RIU,  # Refractive index of the scatterer
        medium_refractive_index=1.0 * ureg.RIU,  # Refractive index of the medium
    )

    detector = CoherentMode(
        mode_number=mode_number,
        NA=0.2 * ureg.AU,  # Numerical aperture for the detector
        sampling=100,  # Field sampling
        gamma_offset=0 * ureg.degree,  # Gamma offset
        phi_offset=0 * ureg.degree,  # Phi offset
        rotation=0 * ureg.degree,  # Rotation angle
    )

    footprint = Footprint(scatterer=scatterer, detector=detector)

    assert footprint is not None, "Expected a valid footprint but got None."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
