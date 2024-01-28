#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure

import PyMieScatt as ps
PyMieScatt_measures = {
    'Qext': 0,
    'Qsca': 1,
    'Qabs': 2,
    'g': 3,
    'Qpr': 4,
    'Qback': 5,
    'Qratio': 6
}


def get_PyMieSim_data(source_set, index, diameters, measure_string: str):
    scatterer_set = Sphere(
        diameter=diameters,
        index=index,
        n_medium=1.,
        source_set=source_set
    )

    experiment = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set,
        detector_set=None
    )

    data = experiment.get(getattr(measure, measure_string))

    return data.y.values.squeeze()


def get_PyMieScatt_data(source_set, index, diameters, measure_string: str):
    PyMieScatt_data = []
    for diameter in diameters:
        efficiencies = ps.MieQ(
            m=index,
            wavelength=source_set.wavelength.values[0],
            diameter=diameter,
        )

        measure_idx = PyMieScatt_measures.get(measure_string)
        data = efficiencies[measure_idx]
        PyMieScatt_data.append(float(data))

    return numpy.asarray(PyMieScatt_data).squeeze()


def get_comparison(wavelength, index, diameters, measure_string: str):
    source_set = Gaussian(
        wavelength=wavelength,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1e-3,
        NA=0.2
    )

    PyMieScatt_data = get_PyMieScatt_data(
        source_set=source_set,
        index=index,
        diameters=diameters,
        measure_string=measure_string
    )

    PyMieSim_data = get_PyMieSim_data(
        source_set=source_set,
        index=index,
        diameters=diameters,
        measure_string=measure_string
    )

    return PyMieSim_data, PyMieScatt_data


@pytest.mark.parametrize('measure_string', ['Qext', 'Qsca', 'Qabs', 'g', 'Qpr'])
def test_comparison(measure_string: str):
    PyMieSim_data, PyMieScatt_data = get_comparison(
        wavelength=632e-9,
        index=1.4 + 0.3j,
        diameters=numpy.geomspace(10e-9, 6000e-9, 800),
        measure_string=measure_string
    )

    discrepency = numpy.isclose(PyMieSim_data, PyMieScatt_data, atol=0, rtol=1e-3)

    if not discrepency.mean() > 0.9:
        raise ValueError('Error: mismatch on PyMieScatt calculation occuring')

# -
