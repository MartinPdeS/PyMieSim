import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian, PolarizationState
from PyMieSim.single.detector import Photodiode


@pytest.mark.parametrize(
    "optical_power, numerical_aperture, polarization, wavelength",
    [
        (1, 0.1 * ureg.AU, 0 * ureg.degree, 1550 * ureg.nanometer),
        (1 * ureg.watt, 0.1 * ureg.AU, 0, 1550 * ureg.nanometer),
        (1 * ureg.watt, 0.1 * ureg.AU, 0 * ureg.degree, 1550),
    ],
)
def test_invalid_gaussian_initialization(optical_power, numerical_aperture, polarization, wavelength):
    with pytest.raises(Exception):
        Gaussian(
            optical_power=optical_power,
            numerical_aperture=numerical_aperture,
            polarization=PolarizationState(angle=polarization),
            wavelength=wavelength,
        )


def test_invalid_sphere_initialization():
    source = Gaussian(
        optical_power=1 * ureg.watt,
        numerical_aperture=0.1 * ureg.AU,
        polarization=PolarizationState(angle=0 * ureg.degree),
        wavelength=1550 * ureg.nanometer,
    )
    with pytest.raises(Exception):
        Sphere(
            source=source,
            medium=1.0 * ureg.RIU,
            material=1.5 * ureg.RIU,
            diameter=100,
        )


@pytest.mark.parametrize(
    "numerical_aperture, sampling, gamma_offset, phi_offset, polarization_filter",
    [
        (0.2 * ureg.AU, 30, 0, 0 * ureg.degree, None),
        (0.2 * ureg.AU, 30, 0 * ureg.degree, 0, None),
    ],
)
def test_invalid_photodiode_initialization(
    numerical_aperture, sampling, gamma_offset, phi_offset, polarization_filter
):
    with pytest.raises(Exception):
        Photodiode(
            numerical_aperture=numerical_aperture,
            sampling=sampling,
            gamma_offset=gamma_offset,
            phi_offset=phi_offset,
            polarization_filter=polarization_filter,
            medium=1.0 * ureg.RIU
        )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
