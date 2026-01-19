#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyMieSim.units import ureg, Angle, Dimensionless
from TypedUnit import AnyUnit
import numpy
from pydantic.dataclasses import dataclass
from typing import List, Union, Optional

from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.experiment.utils import Sequential
from PyMieSim.utils import config_dict
from PyMieSim.binary.interface_experiment import CoherentModeSet

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
    mean_coupling: Optional[bool] = False
    cache_NA: Dimensionless = (0.0,) * ureg.AU
    medium_refractive_index: Dimensionless = (1.0,) * ureg.AU
    sampling: Optional[Dimensionless] = (200,) * ureg.AU
    polarization_filter: Optional[Angle] = (numpy.nan,) * ureg.degree

    def _generate_binding(self) -> None:
        """
        Initializes the C++ binding for the detector using the given simulation parameters. This ensures that the
        detector is correctly linked to the backend, enabling high-performance Mie scattering calculations.

        Sets up parameters such as mode number, sampling rate, NA, and various offsets for the simulation.
        """
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

        self.set = CoherentModeSet(
            **{
                k: v.to_base_units().magnitude if isinstance(v, AnyUnit) else v
                for k, v in self.binding_kwargs.items()
            }
        )