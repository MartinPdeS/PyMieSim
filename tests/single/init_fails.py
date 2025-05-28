import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.units import nanometer, degree, watt, AU, RIU


@pytest.mark.parametrize("optical_power, NA, polarization, wavelength", [
    (1, 0.1 * AU, 0 * degree, 1550 * nanometer),
    (1 * watt, 0.1, 0 * degree, 1550 * nanometer),
    (1 * watt, 0.1 * AU, 0, 1550 * nanometer),
    (1 * watt, 0.1 * AU, 0 * degree, 1550)
])
def test_invalid_gaussian_initialization(optical_power, NA, polarization, wavelength):
    with pytest.raises(ValueError):
        Gaussian(optical_power=optical_power, NA=NA, polarization=polarization, wavelength=wavelength)


@pytest.mark.parametrize("medium_property, property, diameter", [
    (1, 1.5 * RIU, 100 * nanometer),
    (1 * RIU, 1.5, 100 * nanometer),
    (1 * RIU, 1.5 * RIU, 100)
])
def test_invalid_sphere_initialization(medium_property, property, diameter):
    source = Gaussian(optical_power=1 * watt, NA=0.1 * AU, polarization=0 * degree, wavelength=1550 * nanometer)
    with pytest.raises(ValueError):
        Sphere(source=source, medium_property=medium_property, property=property, diameter=diameter)


@pytest.mark.parametrize("NA, sampling, gamma_offset, phi_offset, polarization_filter", [
    (0.2 * AU, 30, 0, 0 * degree, None),
    (0.2 * AU, 30, 0 * degree, 0, None),
    (0.2, 30 * AU, 0 * degree, 0 * degree, None)
])
def test_invalid_photodiode_initialization(NA, sampling, gamma_offset, phi_offset, polarization_filter):
    with pytest.raises(ValueError):
        Photodiode(
            NA=NA,
            sampling=sampling,
            gamma_offset=gamma_offset,
            phi_offset=phi_offset,
            polarization_filter=polarization_filter
        )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
