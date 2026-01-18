#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment import Setup

SIZE = 10
ONES = numpy.ones(SIZE)

WAVELENGTH = numpy.linspace(600, 1000, SIZE) * ureg.nanometer
POLARIZATION = 0 * ureg.degree
OPTICAL_POWER = 1e-3 * ureg.watt
NA = 0.2 * ureg.AU

DIAMETER = 1400 * ureg.nanometer
PROPERTY = 1.4 * ureg.RIU
MEDIUM_PROPERTY = 1.1 * ureg.RIU


CACHE_NA = 0.2 * NA
PHI_OFFSET = 0 * ureg.degree
GAMMA_OFFSET = 0 * ureg.degree
POLARIZATION_FILTER = 0 * ureg.degree
SAMPLING = 100 * ureg.AU
ROTATION = 0 * ureg.degree


@pytest.fixture
def sequential_source():
    return Gaussian(
        wavelength=WAVELENGTH,
        polarization=ONES * POLARIZATION,
        optical_power=ONES * OPTICAL_POWER,
        NA=ONES * NA,
    )


@pytest.fixture
def sequential_scatterer(sequential_source):
    return Sphere(
        diameter=ONES * DIAMETER,
        source=sequential_source,
        property=ONES * PROPERTY,
        medium_property=ONES * MEDIUM_PROPERTY,
    )


@pytest.fixture
def sequential_detector():
    detector = Photodiode(
        phi_offset=ONES * PHI_OFFSET,
        gamma_offset=ONES * GAMMA_OFFSET,
        NA=ONES * NA,
        cache_NA=ONES * CACHE_NA,
        sampling=ONES * SAMPLING,
        polarization_filter=ONES * POLARIZATION_FILTER,
    )

    detector.rotation = ONES * ROTATION

    detector.mode_number = SIZE * ["NC00"]

    detector._generate_binding()

    return detector


@pytest.fixture
def standard_source():
    return Gaussian(
        wavelength=WAVELENGTH,
        polarization=POLARIZATION,
        optical_power=OPTICAL_POWER,
        NA=NA,
    )


@pytest.fixture
def standard_scatterer(standard_source):
    # Configure the spherical scatterer
    return Sphere(
        diameter=DIAMETER,
        source=standard_source,
        property=PROPERTY,
        medium_property=MEDIUM_PROPERTY,
    )


@pytest.fixture
def standard_detector():
    # Configure the spherical scatterer
    return Photodiode(
        phi_offset=PHI_OFFSET,
        gamma_offset=GAMMA_OFFSET,
        NA=NA,
        cache_NA=CACHE_NA,
        sampling=SAMPLING,
        polarization_filter=POLARIZATION_FILTER,
    )


def test_sequential_vs_standard_no_detector(
    sequential_source, sequential_scatterer, standard_source, standard_scatterer
):
    data_standard = Setup(scatterer=standard_scatterer, source=standard_source).get(
        "Qsca", add_units=False
    )

    data_standard = data_standard.values.squeeze()

    assert data_standard is not None, "Error while running the standard get function."

    data_sequential = (
        Setup(scatterer=sequential_scatterer, source=sequential_source)
        .get_sequential("Qsca")
        .squeeze()
    )

    assert data_sequential is not None, "Error while running the standard get function."

    assert not numpy.any(
        data_sequential == 0
    ), "Zeros value retrieved in sequential data means error."

    assert numpy.all(
        data_sequential == data_standard
    ), "Mismatch between sequential and standard computed data."


def test_sequential_vs_standard_detector(
    sequential_source,
    sequential_scatterer,
    sequential_detector,
    standard_source,
    standard_scatterer,
    standard_detector,
):
    data_standard = Setup(
        scatterer=standard_scatterer, source=standard_source, detector=standard_detector
    ).get("coupling", add_units=False)

    data_standard = data_standard.values.squeeze()

    assert data_standard is not None, "Error while running the standard get function."

    data_sequential = (
        Setup(
            scatterer=sequential_scatterer,
            source=sequential_source,
            detector=sequential_detector,
        )
        .get_sequential("coupling")
        .squeeze()
    )

    assert data_sequential is not None, "Error while running the standard get function."

    assert not numpy.any(
        data_sequential == 0
    ), "Zeros value retrieved in sequential data means error."

    assert numpy.all(
        data_sequential == data_standard
    ), "Mismatch between sequential and standard computed data."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
