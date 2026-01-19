#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyMieSim.units import ureg, Angle, Dimensionless
from typing import Optional
import numpy
from dataclasses import field
from pydantic.dataclasses import dataclass
from typing import Tuple
from TypedUnit import AnyUnit
from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.experiment.utils import Sequential
from PyMieSim.utils import config_dict
from PyMieSim.binary.interface_experiment import PhotodiodeSet


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
    cache_NA: Dimensionless = (0.0,) * ureg.AU
    sampling: Optional[Dimensionless] = (200,) * ureg.AU
    polarization_filter: Optional[Angle | None] = (numpy.nan,) * ureg.degree
    mode_number: Tuple[str] = field(default_factory=lambda: ["NC00"], init=True)
    rotation: Angle = field(default=(0,) * ureg.degree, init=True)
    medium_refractive_index: Dimensionless = (1.0,) * ureg.AU
    mean_coupling: bool = field(default=False, init=False)

    def _generate_binding(self) -> None:
        """
        Initializes the C++ binding for the detector using the given simulation parameters. This ensures that the
        detector is correctly linked to the backend, enabling high-performance Mie scattering calculations.

        Sets up parameters such as mode number, sampling rate, NA, and various offsets for the simulation.
        """
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