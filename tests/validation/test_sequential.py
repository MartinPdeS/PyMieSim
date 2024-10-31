#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU

SIZE = 10
ONES = numpy.ones(SIZE)

WAVELENGTH = numpy.linspace(600, 1000, SIZE) * nanometer
POLARIZATION = 0 * degree
OPTICAL_POWER = 1e-3 * watt
NA = 0.2 * AU

DIAMETER = 1400 * nanometer
PROPERTY = 1.4 * RIU
MEDIUM_PROPERTY = 1.1 * RIU


PHI_OFFSET = 0 * degree
GAMMA_OFFSET = 0 * degree
POLARIZATION_FILTER = 0 * degree
SAMPLING = 100 * AU
ROTATION = 0 * degree


@pytest.fixture
def sequential_source():
    return Gaussian(
        wavelength=WAVELENGTH,
        polarization=ONES * POLARIZATION,
        optical_power=ONES * OPTICAL_POWER,
        NA=ONES * NA
    )


@pytest.fixture
def sequential_scatterer(sequential_source):
    return Sphere(
        diameter=ONES * DIAMETER,
        source=sequential_source,
        property=ONES * PROPERTY,
        medium_property=ONES * MEDIUM_PROPERTY
    )


@pytest.fixture
def standard_source():
    return Gaussian(
        wavelength=WAVELENGTH,
        polarization=POLARIZATION,
        optical_power=OPTICAL_POWER,
        NA=NA
    )


@pytest.fixture
def standard_scatterer(standard_source):
    # Configure the spherical scatterer
    return Sphere(
        diameter=DIAMETER,
        source=standard_source,
        property=PROPERTY,
        medium_property=MEDIUM_PROPERTY
    )


def test_sequential_vs_standard_no_detector(sequential_source, sequential_scatterer, standard_source, standard_scatterer):
    data_standard = Setup(
        scatterer=standard_scatterer,
        source=standard_source
    ).get('Qsca', add_units=False)
    print(data_standard)
    data_standard = data_standard.values.squeeze()

    assert data_standard is not None, "Error while running the standard get function."

    data_sequential = Setup(
        scatterer=sequential_scatterer,
        source=sequential_source
    ).get_sequential('Qsca').squeeze()

    assert data_sequential is not None, "Error while running the standard get function."

    assert not numpy.any(data_sequential == 0), "Zeros value retrieved in sequential data means error."

    assert numpy.all(data_sequential == data_standard), "Mismatch between sequential and standard computed data."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
