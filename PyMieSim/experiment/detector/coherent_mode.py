#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import field
from pydantic.dataclasses import dataclass
from typing import List, Union, Optional
from PyMieSim.experiment.detector.base import BaseDetector, config_dict
from PyMieSim.units import Quantity, degree


@dataclass(config=config_dict)
class CoherentMode(BaseDetector):
    """
    Coherent mode detector for Mie scattering simulations, handling coherent detection modes.

    This detector is designed specifically for coherent modes, such as LP, HG, LG, and NC modes, which require specific handling
    in Mie scattering experiments.

    Parameters
    ----------
    mode_number : Union[List[str], str]
        Mode number(s) involved in the detection.
    NA : Union[List[float], float]
        Numerical aperture(s) of the detector.
    gamma_offset : Quantity
        Gamma angular offset (in degrees).
    phi_offset : Quantity
        Phi angular offset (in degrees).
    rotation : Quantity
        Rotation angle of the detector.
    sampling : Union[List[int], int]
        Sampling rate(s) for the detector.
    polarization_filter : Optional[Quantity]
        Polarization filter angle (in degrees).
    mean_coupling : Optional[bool]
        Whether mean coupling is used. Defaults to False.
    coherent : bool
        Specifies if the detection is coherent. Defaults to True.
    """
    rotation: Quantity = 900 * degree
    mode_number: Union[List[str], str]
    mean_coupling: Optional[bool] = False
    coherent: bool = field(default=True, init=False)
