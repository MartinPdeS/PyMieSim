#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment import Setup

SIZE = 10
ONES = numpy.ones(SIZE)

WAVELENGTH = numpy.linspace(600, 1000, SIZE) * ureg.nanometer
POLARIZATION = 0
OPTICAL_POWER = 1e-3 * ureg.watt
NA = 0.2

DIAMETER = 1400 * ureg.nanometer
PROPERTY = 1.4
MEDIUM_PROPERTY = 1.1


CACHE_NA = 0.2 * NA
PHI_OFFSET = 0 * ureg.degree
GAMMA_OFFSET = 0 * ureg.degree
POLARIZATION_FILTER = 0 * ureg.degree
SAMPLING = 100
ROTATION = 0 * ureg.degree


def test_sequential_vs_standard_no_detector():
    polarization_set = PolarizationSet(angles=0 * ureg.degree)
    source_standard = GaussianSet(
        wavelength=WAVELENGTH,
        polarization=polarization_set,
        optical_power=OPTICAL_POWER,
        numerical_aperture=NA,
    )

    scatterer_standard = SphereSet(
        diameter=DIAMETER,
        material=PROPERTY,
        medium=MEDIUM_PROPERTY,
    )

    setup_standard = Setup(
        scatterer_set=scatterer_standard,
        source_set=source_standard
    )

    data_standard = setup_standard.get("Qsca", as_numpy=True).squeeze()

    assert setup_standard is not None, "Error while running the standard get function."

    source_sequential = GaussianSet.build_sequential(
        target_size=SIZE,
        wavelength=WAVELENGTH,
        polarization=PolarizationSet(angles=ONES * ureg.degree),
        optical_power=ONES * OPTICAL_POWER,
        numerical_aperture=ONES * NA,
    )

    scatterer_sequential = SphereSet.build_sequential(
        target_size=SIZE,
        diameter=ONES * DIAMETER,
        material=ONES * PROPERTY,
        medium=ONES * MEDIUM_PROPERTY,
    )

    setup_sequential = Setup(
        scatterer_set=scatterer_sequential,
        source_set=source_sequential
    )

    data_sequential = setup_sequential.get_sequential("Qsca").squeeze()

    assert data_sequential is not None, "Error while running the standard get function."

    assert not numpy.any(
        data_sequential == 0
    ), "Zeros value retrieved in sequential data means error."

    assert numpy.all(
        data_sequential == data_standard
    ), "Mismatch between sequential and standard computed data."


def test_sequential_vs_standard_detector():
    source_standard = GaussianSet(
        wavelength=WAVELENGTH,
        polarization=0 * ureg.degree,
        optical_power=OPTICAL_POWER,
        numerical_aperture=NA,
    )

    scatterer_standard = SphereSet(
        diameter=DIAMETER,
        material=PROPERTY,
        medium=MEDIUM_PROPERTY,
    )

    detector_standard = PhotodiodeSet(
        phi_offset=PHI_OFFSET,
        gamma_offset=GAMMA_OFFSET,
        numerical_aperture=NA,
        cache_numerical_aperture=CACHE_NA,
        sampling=SAMPLING,
        polarization_filter=0 * ureg.degree
    )

    setup_standard = Setup(
        scatterer_set=scatterer_standard,
        source_set=source_standard,
        detector_set=detector_standard
    )

    data_standard = setup_standard.get("coupling", as_numpy=True).squeeze()

    assert data_standard is not None, "Error while running the standard get function."


    source_sequential = GaussianSet.build_sequential(
        target_size=SIZE,
        wavelength=WAVELENGTH,
        polarization=PolarizationSet(angles=POLARIZATION * ONES * ureg.degree),
        optical_power=ONES * OPTICAL_POWER,
        numerical_aperture=ONES * NA,
    )

    scatterer_sequential = SphereSet.build_sequential(
        target_size=SIZE,
        diameter=ONES * DIAMETER,
        material=ONES * PROPERTY,
        medium=ONES * MEDIUM_PROPERTY,
    )

    detector_sequential = PhotodiodeSet.build_sequential(
        target_size=SIZE,
        phi_offset=ONES * PHI_OFFSET,
        gamma_offset=ONES * GAMMA_OFFSET,
        numerical_aperture=ONES * NA,
        cache_numerical_aperture=ONES * CACHE_NA,
        sampling=[SAMPLING] * SIZE,
        polarization_filter=ONES * POLARIZATION_FILTER,
        medium=ONES * MEDIUM_PROPERTY,
    )

    setup_sequential = Setup(
        scatterer_set=scatterer_sequential,
        source_set=source_sequential,
        detector_set=detector_sequential,
    )

    data_sequential = setup_sequential.get_sequential("coupling")

    assert data_sequential is not None, "Error while running the standard get function."

    assert not numpy.any(
        data_sequential == 0
    ), "Zeros value retrieved in sequential data means error."

    assert numpy.all(
        data_sequential == data_standard
    ), "Mismatch between sequential and standard computed data."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
