import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian, PolarizationState
from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single.representations import Footprint

# Define sampling to be tested
samplings = [100, 200]


@pytest.mark.parametrize("sampling", samplings)
def test_photodiode_with_sampling(sampling):
    """Test the Integrating Sphere detector with different sampling rates."""

    source = Gaussian(
        wavelength=750 * ureg.nanometer,  # Wavelength of the source in meters
        polarization=PolarizationState(angle=0 * ureg.degree),  # Polarization value
        optical_power=1 * ureg.watt,  # Optical power in watts
        numerical_aperture=0.3 * ureg.AU,  # Numerical aperture
    )

    scatterer =  Sphere(
        diameter=100 * ureg.nanometer,  # Diameter of the scatterer in meters
        source=source,  # Gaussian source from source fixture
        material=1.4 * ureg.RIU,  # Refractive index of the scatterer
        medium=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    )

    detector = IntegratingSphere(sampling=sampling)

    footprint = Footprint(scatterer=scatterer, detector=detector)

    # Example verification step (not operational as we're not evaluating the output here)
    assert footprint is not None, "Expected a valid footprint but got None."

if __name__ == "__main__":
    pytest.main(["-W error", __file__])
