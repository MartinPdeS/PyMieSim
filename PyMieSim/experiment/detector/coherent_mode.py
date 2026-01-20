#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from typing import List, Union, Optional

from PyMieSim.units import ureg, Angle, Dimensionless
from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.experiment.utils import Sequential
from PyMieSim.binary.interface_experiment import CoherentModeSet


class CoherentMode(BaseDetector, Sequential):
    """
    Coherent mode detector for Mie scattering simulations, handling coherent detection modes.

    This detector is designed specifically for coherent modes, such as LP, HG, LG, and NC modes, which require specific handling
    in Mie scattering experiments.

    Parameters
    ----------
    mode_number : Union[List[str], str]
        Mode number(s) involved in the detection.
    NA : Dimensionless
        Numerical aperture(s) of the detector.
    gamma_offset : Angle
        Gamma angular offset (in degrees).
    phi_offset : Angle
        Phi angular offset (in degrees).
    rotation : Angle
        Rotation angle of the detector.
    sampling : Dimensionless
        Sampling rate(s) for the detector.
    polarization_filter : Optional[Angle]
        Polarization filter angle (in degrees).
    mean_coupling : Optional[bool]
        Whether mean coupling is used. Defaults to False.
    coherent : bool
        Specifies if the detection is coherent. Defaults to True.
    """

    attributes = [
        "mode_number",
        "NA",
        "cache_NA",
        "gamma_offset",
        "phi_offset",
        "rotation",
        "sampling",
        "polarization_filter",
        "medium_refractive_index",
    ]

    def __init__(
        self,
        mode_number: Union[List[str], str],
        NA: Dimensionless,
        gamma_offset: Angle,
        phi_offset: Angle,
        rotation: Angle,
        mean_coupling: Optional[bool] = False,
        cache_NA: Dimensionless = (0.0,) * ureg.AU,
        medium_refractive_index: Dimensionless = (1.0,) * ureg.AU,
        sampling: Optional[Dimensionless] = (200,) * ureg.AU,
        polarization_filter: Optional[Angle] = (numpy.nan,) * ureg.degree
    ):
        self.mode_number = numpy.atleast_1d(mode_number)

        for mode in self.mode_number:
            assert mode[:2] in ['LP', 'HG', 'LG', 'NC'], f"mode_number must be one of 'LP', 'HG', 'LG', or 'NC', got {mode_number}"
        self.NA = numpy.atleast_1d(NA)
        self.gamma_offset = numpy.atleast_1d(gamma_offset)
        self.phi_offset = numpy.atleast_1d(phi_offset)
        self.rotation = numpy.atleast_1d(rotation)
        self.sampling = numpy.atleast_1d(sampling)
        self.cache_NA = numpy.atleast_1d(cache_NA)
        self.medium_refractive_index = numpy.atleast_1d(medium_refractive_index)
        self.mean_coupling = mean_coupling


        if polarization_filter is None:
            polarization_filter = numpy.nan * ureg.degree

        self.polarization_filter = numpy.atleast_1d(polarization_filter).astype(float)

        self.binding_kwargs = {
            "mode_number": self.mode_number,
            "sampling": self.sampling,
            "NA": self.NA,
            "cache_NA": self.cache_NA,
            "phi_offset": self.phi_offset,
            "gamma_offset": self.gamma_offset,
            "polarization_filter": self.polarization_filter,
            "rotation": self.rotation,
            "medium_refractive_index": self.medium_refractive_index,
            "mean_coupling": self.mean_coupling,
            "is_sequential": self.is_sequential,
        }

        self.set = CoherentModeSet(**self.binding_kwargs)