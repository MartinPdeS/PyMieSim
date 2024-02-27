#!/usr/bin/env python
# -*- coding: utf-8 -*-


from PyMieSim import single
from PyMieSim import experiment, measure


def test_detector_single_polarization_filter():
    source = single.source.Gaussian(
        wavelength=1e-6,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )

    scatterer = single.scatterer.Sphere(
        diameter=1e-6,
        index=1.5 + 0.5j,
        source=source
    )

    detector_0 = single.detector.Photodiode(
        NA=0.1,
        gamma_offset=0,
        phi_offset=90,
        polarization_filter=0
    )

    detector_180 = single.detector.Photodiode(
        NA=0.1,
        gamma_offset=0,
        phi_offset=90,
        polarization_filter=180
    )

    if not detector_0.coupling(scatterer) == detector_180.coupling(scatterer):
        raise ValueError('Mismatch with coupling value for detector with polarization filter of 0 and 180 degrees')


def test_detector_experiment_polarization_filter():
    source_set = experiment.source.Gaussian(
        wavelength=1e-6,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )

    scatterer = experiment.scatterer.Sphere(
        n_medium=1.0,
        diameter=1e-6,
        index=1.5 + 0.5j,
        source_set=source_set
    )

    detector = experiment.detector.Photodiode(
        NA=0.1,
        gamma_offset=0,
        phi_offset=90,
        polarization_filter=[0, 180],
        sampling=100
    )

    setup = experiment.Setup(
        scatterer_set=scatterer,
        detector_set=detector,
        source_set=source_set
    )

    coupling_values = setup.get(measure=measure.coupling, export_as_numpy=True).squeeze()

    if not coupling_values[0] == coupling_values[-1]:
        raise ValueError('Mismatch with coupling value for detector with polarization filter of 0 and 180 degrees')

# -
