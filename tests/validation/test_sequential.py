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




def test_sequential_vs_standard_no_detector():
    source_standard = Gaussian(
        wavelength=WAVELENGTH,
        polarization=POLARIZATION,
        optical_power=OPTICAL_POWER,
        NA=NA,
    )

    scatterer_standard = Sphere(
        diameter=DIAMETER,
        source=source_standard,
        refractive_index=PROPERTY,
        medium_refractive_index=MEDIUM_PROPERTY,
    )

    setup_standard = Setup(
        scatterer=scatterer_standard,
        source=source_standard
    )

    data_standard = setup_standard.get("Qsca", add_units=False).values.squeeze()

    assert setup_standard is not None, "Error while running the standard get function."

    source_sequential = Gaussian(
        wavelength=WAVELENGTH,
        polarization=ONES * POLARIZATION,
        optical_power=ONES * OPTICAL_POWER,
        NA=ONES * NA,
    )

    scatterer_sequential = Sphere(
        diameter=ONES * DIAMETER,
        source=source_sequential,
        refractive_index=ONES * PROPERTY,
        medium_refractive_index=ONES * MEDIUM_PROPERTY,
    )

    setup_sequential = Setup(
        scatterer=scatterer_sequential,
        source=source_sequential
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
    source_standard = Gaussian(
        wavelength=WAVELENGTH,
        polarization=POLARIZATION,
        optical_power=OPTICAL_POWER,
        NA=NA,
    )

    scatterer_standard = Sphere(
        diameter=DIAMETER,
        source=source_standard,
        refractive_index=PROPERTY,
        medium_refractive_index=MEDIUM_PROPERTY,
    )

    detector_standard = Photodiode(
        phi_offset=PHI_OFFSET,
        gamma_offset=GAMMA_OFFSET,
        NA=NA,
        cache_NA=CACHE_NA,
        sampling=SAMPLING,
        polarization_filter=POLARIZATION_FILTER,
    )

    setup_standard = Setup(
        scatterer=scatterer_standard,
        source=source_standard,
        detector=detector_standard
    )

    data_standard = setup_standard.get("coupling", add_units=False).values.squeeze()

    assert data_standard is not None, "Error while running the standard get function."


    source_sequential = Gaussian(
        wavelength=WAVELENGTH,
        polarization=ONES * POLARIZATION,
        optical_power=ONES * OPTICAL_POWER,
        NA=ONES * NA,
    )

    scatterer_sequential = Sphere(
        diameter=ONES * DIAMETER,
        source=source_sequential,
        refractive_index=ONES * PROPERTY,
        medium_refractive_index=ONES * MEDIUM_PROPERTY,
    )

    detector_sequential = Photodiode(
        phi_offset=ONES * PHI_OFFSET,
        gamma_offset=ONES * GAMMA_OFFSET,
        NA=ONES * NA,
        cache_NA=ONES * CACHE_NA,
        sampling=ONES * SAMPLING,
        polarization_filter=ONES * POLARIZATION_FILTER,
        medium_refractive_index=ONES * MEDIUM_PROPERTY,
    )

    setup_sequential = Setup(
        scatterer=scatterer_sequential,
        source=source_sequential,
        detector=detector_sequential,
    )

    data_sequential = setup_sequential.get_sequential("coupling").squeeze()

    assert data_sequential is not None, "Error while running the standard get function."

    assert not numpy.any(
        data_sequential == 0
    ), "Zeros value retrieved in sequential data means error."

    assert numpy.all(
        data_sequential == data_standard
    ), "Mismatch between sequential and standard computed data."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
