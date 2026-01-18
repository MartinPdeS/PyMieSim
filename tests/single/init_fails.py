import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode


@pytest.mark.parametrize(
    "optical_power, NA, polarization, wavelength",
    [
        (1, 0.1 * ureg.AU, 0 * ureg.degree, 1550 * ureg.nanometer),
        (1 * ureg.watt, 0.1 * ureg.AU, 0, 1550 * ureg.nanometer),
        (1 * ureg.watt, 0.1 * ureg.AU, 0 * ureg.degree, 1550),
    ],
)
def test_invalid_gaussian_initialization(optical_power, NA, polarization, wavelength):
    with pytest.raises(Exception):
        Gaussian(
            optical_power=optical_power,
            NA=NA,
            polarization=polarization,
            wavelength=wavelength,
        )


@pytest.mark.parametrize(
    "medium_property, property, diameter",
    [
        (1, 1.5 * ureg.RIU, 100 * ureg.nanometer),
        (1 * ureg.RIU, 1.5, 100 * ureg.nanometer),
        (1 * ureg.RIU, 1.5 * ureg.RIU, 100),
    ],
)
def test_invalid_sphere_initialization(medium_property, property, diameter):
    source = Gaussian(
        optical_power=1 * ureg.watt,
        NA=0.1 * ureg.AU,
        polarization=0 * ureg.degree,
        wavelength=1550 * ureg.nanometer,
    )
    with pytest.raises(Exception):
        Sphere(
            source=source,
            medium_property=medium_property,
            property=property,
            diameter=diameter,
        )


@pytest.mark.parametrize(
    "NA, sampling, gamma_offset, phi_offset, polarization_filter",
    [
        (0.2 * ureg.AU, 30, 0, 0 * ureg.degree, None),
        (0.2 * ureg.AU, 30, 0 * ureg.degree, 0, None),
    ],
)
def test_invalid_photodiode_initialization(
    NA, sampling, gamma_offset, phi_offset, polarization_filter
):
    with pytest.raises(Exception):
        Photodiode(
            NA=NA,
            sampling=sampling,
            gamma_offset=gamma_offset,
            phi_offset=phi_offset,
            polarization_filter=polarization_filter,
            medium_refractive_index=1.0 * ureg.RIU
        )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
