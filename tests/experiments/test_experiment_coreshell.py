#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy
import PyMieSim
from PyMieSim.experiment import CoreShellSet, SourceSet, PhotodiodeSet, Setup
from PyMieSim.materials import Silver, BK7, Aluminium

core_type = [
    {'name': 'BK7', 'kwarg': {'core_material': BK7}},
    {'name': 'Silver', 'kwarg': {'core_material': Silver}},
    {'name': 'Aluminium', 'kwarg': {'core_material': Aluminium}},
    {'name': 'Index', 'kwarg': {'core_index': 1.4}}
]

shell_type = [
    {'name': 'BK7', 'kwarg': {'shell_material': BK7}},
    {'name': 'Silver', 'kwarg': {'shell_material': Silver}},
    {'name': 'Aluminium', 'kwarg': {'shell_material': Aluminium}},
    {'name': 'Index', 'kwarg': {'shell_index': 1.4}}
]

measures = [
    PyMieSim.measure.Qsca,
    PyMieSim.measure.Qabs,
    PyMieSim.measure.Qback,
    PyMieSim.measure.g,
    PyMieSim.measure.a1,
    PyMieSim.measure.b1,
    PyMieSim.measure.coupling
]


@pytest.mark.parametrize('shell_type', [p['kwarg'] for p in core_type], ids=[p['name'] for p in core_type])
@pytest.mark.parametrize('core_type', [p['kwarg'] for p in shell_type], ids=[p['name'] for p in shell_type])
@pytest.mark.parametrize('measure', measures, ids=[p.name for p in measures])
def test_coreshell_experiment(measure, core_type, shell_type):
    scatterer_set = CoreShellSet(
        n_medium=1,
        core_diameter=numpy.linspace(400e-9, 1400e-9, 10),
        shell_width=300e-9,
        **core_type,
        **shell_type,
    )

    source_set = SourceSet(
        wavelength=numpy.linspace(400e-9, 1800e-9, 50),
        polarization=0,
        amplitude=1
    )

    detector_set = PhotodiodeSet(
        NA=0.2,
        polarization_filter=None,
        gamma_offset=0,
        phi_offset=0
    )

    experiment = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set,
        detector_set=detector_set
    )

    experiment.Get(measure)

# -
