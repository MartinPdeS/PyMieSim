#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from dataclasses import field
from pydantic.dataclasses import dataclass
from pydantic import validator
from PyMieSim.units import Quantity
from typing import List, Union, Optional
from PyMieSim.experiment.detector.base import BaseDetector, config_dict

@dataclass(config=config_dict)
class CoherentMode(BaseDetector):
    """
    Coherent mode detector for Mie scattering simulations, handling coherent detection modes.

    This detector is designed specifically for coherent modes, such as LP, HG, LG, and NC modes, which require specific handling
    in Mie scattering experiments.

    Attributes:
        mode_number (Union[List[str], str]): Mode number(s) involved in the detection.
        NA (Union[List[float], float]): Numerical aperture(s) of the detector.
        gamma_offset (Quantity): Gamma angular offset (in degrees).
        phi_offset (Quantity): Phi angular offset (in degrees).
        rotation (Quantity): Rotation angle of the detector.
        sampling (Union[List[int], int]): Sampling rate(s) for the detector.
        polarization_filter (Optional[Quantity]): Polarization filter angle (in degrees).
        mean_coupling (Optional[bool]): Whether mean coupling is used. Defaults to False.
        coherent (bool): Specifies if the detection is coherent. Defaults to True.
    """
    mode_number: Union[np.ndarray, List[str], str]
    NA: Quantity
    gamma_offset: Quantity
    phi_offset: Quantity
    rotation: Quantity
    sampling: Quantity
    polarization_filter: Optional[Quantity] = None
    mean_coupling: Optional[bool] = False

    coherent: bool = field(default=True, init=False)

    @validator('mode_number', pre=True)
    def validate_mode_number(cls, mode_number):
        """Ensure mode numbers are valid and belong to supported families."""
        mode_number = np.atleast_1d(mode_number).astype(str)
        for mode in mode_number:
            if mode[:2] not in ['LP', 'HG', 'LG', 'NC']:
                raise ValueError(f'Invalid mode family {mode[:2]}. Must be one of: LP, HG, LG, NC')
        return mode_number
