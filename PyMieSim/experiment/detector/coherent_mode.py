#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from typing import List, Union, Optional
from TypedUnit import ureg, Angle, Dimensionless

from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.experiment.utils import Sequential
from PyMieSim.utils import config_dict


@dataclass(config=config_dict)
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

    mode_number: Union[List[str], str]
    NA: Dimensionless
    gamma_offset: Angle
    phi_offset: Angle
    rotation: Angle
    mean_coupling: bool
    coherent: bool = True
    mean_coupling: Optional[bool] = False
    cache_NA: Dimensionless = (0.0,) * ureg.AU
    sampling: Optional[Dimensionless] = (200,) * ureg.AU
    polarization_filter: Optional[Angle] = (numpy.nan,) * ureg.degree
