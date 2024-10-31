#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU, volt, meter


@pytest.mark.parametrize("wavelength, polarization, optical_power, NA", [
    (100, 0 * degree, 1e-3 * watt, 0.2 * AU),
    (100 * nanometer, 0 * degree, 1e-3, 0.2 * AU),
    (100 * nanometer, 0 * degree, 1e-3 * watt, 0.2),
    (100 * nanometer, 0, 1e-3 * watt, 0.2 * AU)
])
def test_invalid_gaussian_initialization(wavelength, polarization, optical_power, NA):
    with pytest.raises(ValueError):
        Gaussian(wavelength=wavelength, polarization=polarization, optical_power=optical_power, NA=NA)


@pytest.mark.parametrize("diameter, property, medium_property", [
    (100, 1.5 * RIU, Material.water),
    (100 * nanometer, 1.5, Material.water),
    (100 * nanometer, 1.5 * RIU, 1.0)
])
def test_invalid_sphere_initialization(diameter, property, medium_property):
    source = PlaneWave(wavelength=1e3 * nanometer, polarization=0 * degree, amplitude=1 * volt / meter)
    with pytest.raises(ValueError):
        Sphere(diameter=diameter, source=source, property=property, medium_property=medium_property)


@pytest.mark.parametrize("mode_number, rotation, NA, polarization_filter, gamma_offset, phi_offset, sampling", [
    ('invalid', 0 * degree, 0.2 * AU, None, 0 * degree, 0 * degree, 100 * AU),
    ('LP01', 0, 0.2 * AU, None, 0 * degree, 0 * degree, 100 * AU),
    ('LP01', 0 * degree, 0.2, None, 0 * degree, 0 * degree, 100 * AU),
    ('LP01', 0 * degree, 0.2 * AU, 10, 0 * degree, 0 * degree, 100 * AU),
    ('LP01', 0 * degree, 0.2 * AU, None, 0, 0 * degree, 100 * AU),
    ('LP01', 0 * degree, 0.2 * AU, None, 0 * degree, 0, 100 * AU),
    ('LP01', 0 * degree, 0.2 * AU, None, 0 * degree, 0 * degree, 100)
])
def test_invalid_coherent_mode_initialization(mode_number, rotation, NA, polarization_filter, gamma_offset, phi_offset, sampling):
    with pytest.raises(ValueError):
        CoherentMode(
            mode_number=mode_number,
            rotation=rotation,
            NA=NA,
            polarization_filter=polarization_filter,
            gamma_offset=gamma_offset,
            phi_offset=phi_offset,
            sampling=sampling
        )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
