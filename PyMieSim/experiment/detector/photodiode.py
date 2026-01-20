#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Optional
import numpy

from PyMieSim.units import ureg, Angle, Dimensionless
from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.experiment.utils import Sequential
from PyMieSim.binary.interface_experiment import PhotodiodeSet


class Photodiode(BaseDetector, Sequential):
    """
    Photodiode detector tailored for Mie scattering simulations, extending BaseDetector with specific features.

    Parameters
    ----------
    NA : Dimensionless
        Numerical aperture(s) of the detector.
    gamma_offset : Angle
        Gamma angular offset (in degrees).
    phi_offset : Angle
        Phi angular offset (in degrees).
    polarization_filter : Optional[Angle]
        Polarization filter angle (in degrees).
    sampling : Dimensionless
        Sampling rate(s) for the detector.
    coherent : bool
        Indicates if the detection is coherent. Defaults to False.
    """

    attributes = [
        "NA",
        "gamma_offset",
        "phi_offset",
        "sampling",
        "polarization_filter",
        "medium_refractive_index",
    ]

    def __init__(
        self,
        NA: Dimensionless,
        gamma_offset: Angle,
        phi_offset: Angle,
        cache_NA: Dimensionless = (0.0,) * ureg.AU,
        sampling: Optional[Dimensionless] = (200,) * ureg.AU,
        polarization_filter: Optional[Angle | None] = (numpy.nan,) * ureg.degree,
        medium_refractive_index: Dimensionless = (1.0,) * ureg.AU,
    ):
        self.NA = numpy.atleast_1d(NA)
        self.gamma_offset = numpy.atleast_1d(gamma_offset)
        self.phi_offset = numpy.atleast_1d(phi_offset)
        self.cache_NA = numpy.atleast_1d(cache_NA)
        self.medium_refractive_index = numpy.atleast_1d(medium_refractive_index)
        self.sampling = numpy.atleast_1d(sampling)

        if polarization_filter is None:
            polarization_filter = numpy.nan * ureg.degree

        self.polarization_filter = numpy.atleast_1d(polarization_filter).astype(float)


        self.binding_kwargs = {
            "sampling": self.sampling,
            "NA": self.NA,
            "cache_NA": self.cache_NA,
            "phi_offset": self.phi_offset,
            "gamma_offset": self.gamma_offset,
            "polarization_filter": self.polarization_filter,
            "medium_refractive_index": self.medium_refractive_index,
            "is_sequential": self.is_sequential,
        }

        self.set = PhotodiodeSet(**self.binding_kwargs)