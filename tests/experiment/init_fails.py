#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg
from PyOptik import Material

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian, PlaneWave


@pytest.mark.parametrize(
    "wavelength, polarization, optical_power, NA",
    [
        (100, 0 * ureg.degree, 1e-3 * ureg.watt, 0.2 * ureg.AU),
        (100 * ureg.nanometer, 0 * ureg.degree, 1e-3, 0.2 * ureg.AU),
        (100 * ureg.nanometer, 0 * ureg.degree, 1e-3 * ureg.watt, 0.2),
        (100 * ureg.nanometer, 0, 1e-3 * ureg.watt, 0.2 * ureg.AU),
    ],
)
def test_invalid_gaussian_initialization(wavelength, polarization, optical_power, NA):
    with pytest.raises(ValueError):
        Gaussian(
            wavelength=wavelength,
            polarization=polarization,
            optical_power=optical_power,
            NA=NA,
        )


@pytest.mark.parametrize(
    "diameter, refractive_index, medium_refractive_index",
    [
        (100, 1.5 * ureg.RIU, Material.water),
        (100 * ureg.nanometer, 1.5, Material.water),
        (100 * ureg.nanometer, 1.5 * ureg.RIU, 1.0),
    ],
)
def test_invalid_sphere_initialization(diameter, refractive_index, medium_refractive_index):
    source = PlaneWave(
        wavelength=1e3 * ureg.nanometer,
        polarization=0 * ureg.degree,
        amplitude=1 * ureg.volt / ureg.meter,
    )
    with pytest.raises(ValueError):
        Sphere(
            diameter=diameter,
            source=source,
            refractive_index=refractive_index,
            medium_refractive_index=medium_refractive_index,
        )


@pytest.mark.parametrize(
    "mode_number, rotation, NA, polarization_filter, gamma_offset, phi_offset, sampling",
    [
        (
            "invalid",
            0 * ureg.degree,
            0.2 * ureg.AU,
            None,
            0 * ureg.degree,
            0 * ureg.degree,
            100 * ureg.AU,
        ),
        (
            "LP01",
            0,
            0.2 * ureg.AU,
            None,
            0 * ureg.degree,
            0 * ureg.degree,
            100 * ureg.AU,
        ),
        (
            "LP01",
            0 * ureg.degree,
            0.2,
            None,
            0 * ureg.degree,
            0 * ureg.degree,
            100 * ureg.AU,
        ),
        (
            "LP01",
            0 * ureg.degree,
            0.2 * ureg.AU,
            10,
            0 * ureg.degree,
            0 * ureg.degree,
            100 * ureg.AU,
        ),
        (
            "LP01",
            0 * ureg.degree,
            0.2 * ureg.AU,
            None,
            0,
            0 * ureg.degree,
            100 * ureg.AU,
        ),
        (
            "LP01",
            0 * ureg.degree,
            0.2 * ureg.AU,
            None,
            0 * ureg.degree,
            0,
            100 * ureg.AU,
        ),
        (
            "LP01",
            0 * ureg.degree,
            0.2 * ureg.AU,
            None,
            0 * ureg.degree,
            0 * ureg.degree,
            100,
        ),
    ],
)
def test_invalid_coherent_mode_initialization(
    mode_number, rotation, NA, polarization_filter, gamma_offset, phi_offset, sampling
):
    with pytest.raises(ValueError):
        CoherentMode(
            mode_number=mode_number,
            rotation=rotation,
            NA=NA,
            polarization_filter=polarization_filter,
            gamma_offset=gamma_offset,
            phi_offset=phi_offset,
            sampling=sampling,
        )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
