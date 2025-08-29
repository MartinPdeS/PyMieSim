#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Optional
import numpy
from dataclasses import field
from pydantic.dataclasses import dataclass
from typing import Tuple
from TypedUnit import ureg, Angle, Dimensionless

from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.experiment.utils import Sequential
from PyMieSim.utils import config_dict


@dataclass(config=config_dict)
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
    mean_coupling : bool
        Whether mean coupling is used. Defaults to True.
    rotation : Quantity
        Rotation angle of the detector. Defaults to 0 degrees.
    coherent : bool
        Indicates if the detection is coherent. Defaults to False.
    mode_number : str
        Mode number of the detector. Defaults to 'NC00'.
    """

    NA: Dimensionless
    gamma_offset: Angle
    phi_offset: Angle
    mean_coupling: bool
    coherent: bool
    cache_NA: Dimensionless = (0.0,) * ureg.AU
    sampling: Optional[Dimensionless] = (200,) * ureg.AU
    polarization_filter: Optional[Angle | None] = (numpy.nan,) * ureg.degree
    mode_number: Tuple[str] = field(default_factory=lambda: ["NC00"], init=True)
    rotation: Angle = field(default=(0,) * ureg.degree, init=True)

    coherent: bool = field(default=False, init=False)
    mean_coupling: bool = field(default=False, init=False)
