#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import field
from pydantic.dataclasses import dataclass
from PyMieSim.units import Quantity, degree
from typing import Tuple
from PyMieSim.experiment.detector.base import BaseDetector, config_dict


@dataclass(config=config_dict)
class Photodiode(BaseDetector):
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
    coherent: bool = field(default=False, init=False)
    mean_coupling: bool = field(default=False, init=False)
    mode_number: Tuple[str] = field(default_factory=lambda: ['NC00'], init=False)
    rotation: Quantity = field(default=(0,) * degree, init=False)
