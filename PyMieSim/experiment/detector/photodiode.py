#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Optional
import numpy
from dataclasses import field
from pydantic.dataclasses import dataclass
from PyMieSim.units import Quantity, degree, AU
from typing import Tuple
from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.experiment.utils import config_dict, Sequential


@dataclass(config=config_dict)
class Photodiode(BaseDetector, Sequential):
    """
    Photodiode detector tailored for Mie scattering simulations, extending BaseDetector with specific features.

    Parameters
    ----------
    NA : Union[List[float], float]
        Numerical aperture(s) of the detector.
    gamma_offset : Quantity
        Gamma angular offset (in degrees).
    phi_offset : Quantity
        Phi angular offset (in degrees).
    polarization_filter : Optional[Quantity]
        Polarization filter angle (in degrees).
    sampling : Union[List[int], int]
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
    NA: Quantity
    gamma_offset: Quantity
    phi_offset: Quantity
    mean_coupling: bool
    coherent: bool
    cache_NA: Quantity = (0.,) * AU
    sampling: Optional[Quantity] = (200,) * AU
    polarization_filter: Optional[Quantity | None] = (numpy.nan, ) * degree
    mode_number: Tuple[str] = field(default_factory=lambda: ['NC00'], init=True)
    rotation: Quantity = field(default=(0,) * degree, init=True)

    coherent: bool = field(default=False, init=False)
    mean_coupling: bool = field(default=False, init=False)

