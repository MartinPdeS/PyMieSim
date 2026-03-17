#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.units import ureg
from PyOptik import Material

from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet, PlaneWaveSet
from PyMieSim.experiment.polarization_set import PolarizationSet



def test_gaussian_rejects_invalid_wavelength_type():
    with pytest.raises(ValueError):
        GaussianSet(
            wavelength=100,  # must carry units
            polarization=PolarizationSet(angles=0 * ureg.degree),
            optical_power=1e-3 * ureg.watt,
            numerical_aperture=0.2,
        )


def test_gaussian_rejects_invalid_optical_power_units():
    with pytest.raises(ValueError):
        GaussianSet(
            wavelength=100 * ureg.nanometer,
            polarization=PolarizationSet(angles=0 * ureg.degree),
            optical_power=1e-3,  # must carry watt
            numerical_aperture=0.2,
        )

def test_gaussian_rejects_invalid_polarization_type():
    with pytest.raises(TypeError):
        GaussianSet(
            wavelength=100 * ureg.nanometer,
            polarization=0,  # must carry degree
            optical_power=1e-3 * ureg.watt,
            numerical_aperture=0.2,
        )


def _valid_plane_wave_source():
    return PlaneWaveSet(
        wavelength=1e3 * ureg.nanometer,
        polarization=PolarizationSet(angles=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    )


def test_sphere_rejects_invalid_diameter_type():
    with pytest.raises((ValueError, RuntimeError)):
        SphereSet(
            diameter=100,  # must carry length units
            material=1.5 * ureg.RIU,
            medium=Material.water,
        )


def test_sphere_rejects_invalid_refractive_index_type():
    with pytest.raises((AttributeError, RuntimeError)):
        SphereSet(
            diameter=100 * ureg.nanometer,
            material=1.5,  # must carry RIU
            medium=Material.water,
        )


def test_sphere_rejects_invalid_medium_refractive_index_type():
    with pytest.raises((RuntimeError, AttributeError, TypeError)):
        SphereSet(
            diameter=100 * ureg.nanometer,
            material=1.5 * ureg.RIU,
            medium=1.0,  # must be RIU or a Material (per your API)
        )


def test_coherent_mode_rejects_invalid_mode_string():
    with pytest.raises(Exception):
        CoherentModeSet(
            mode_number="invalid",  # should be like "LP01"
            rotation=0 * ureg.degree,
            numerical_aperture=0.2,
            polarization_filter=None,
            gamma_offset=0 * ureg.degree,
            phi_offset=0 * ureg.degree,
            sampling=100,
        )


def test_coherent_mode_rejects_rotation_without_units():
    with pytest.raises(ValueError):
        CoherentModeSet(
            mode_number="LP01",
            rotation=0,  # must carry degree
            numerical_aperture=0.2,
            polarization_filter=None,
            gamma_offset=0 * ureg.degree,
            phi_offset=0 * ureg.degree,
            sampling=100,
        )


def test_coherent_mode_rejects_polarization_filter_wrong_type():
    with pytest.raises(ValueError):
        CoherentModeSet(
            mode_number="LP01",
            rotation=0 * ureg.degree,
            numerical_aperture=0.2,
            polarization_filter=10,  # expected None or an angle with units
            gamma_offset=0 * ureg.degree,
            phi_offset=0 * ureg.degree,
            sampling=100,
        )


def test_coherent_mode_rejects_gamma_offset_without_units():
    with pytest.raises(ValueError):
        CoherentModeSet(
            mode_number="LP01",
            rotation=0 * ureg.degree,
            numerical_aperture=0.2,
            polarization_filter=None,
            gamma_offset=0,  # must carry degree
            phi_offset=0 * ureg.degree,
            sampling=100,
        )


def test_coherent_mode_rejects_phi_offset_without_units():
    with pytest.raises(ValueError):
        CoherentModeSet(
            mode_number="LP01",
            rotation=0 * ureg.degree,
            numerical_aperture=0.2,
            polarization_filter=None,
            gamma_offset=0 * ureg.degree,
            phi_offset=0,  # must carry degree
            sampling=100,
        )


if __name__ == "__main__":
    pytest.main(["-W", "error", __file__])
