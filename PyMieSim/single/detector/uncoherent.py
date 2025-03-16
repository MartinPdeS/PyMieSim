#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from typing import Optional
from dataclasses import field
from pydantic.dataclasses import dataclass
from PyMieSim.binary.interface_detector import BindedDetector
from PyMieSim.units import Quantity, degree, AU, radian
from PyMieSim.single.detector.base import BaseDetector, config_dict


@dataclass(config=config_dict)
class Photodiode(BaseDetector):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This means it is independent of the phase of the impinging scattered light field.

    Parameters
    ----------
    NA : Quantity
        Numerical aperture of the imaging system.
    cache_NA : Quantity
        Numerical aperture of the detector cache.
    gamma_offset : Quantity
        Angle [Degree] offset of the detector in the direction perpendicular to polarization.
    phi_offset : Quantity
        Angle [Degree] offset of the detector in the direction parallel to polarization.
    sampling : Quantity
        Sampling rate of the far-field distribution. Default is 200.
    polarization_filter : Optional[Quantity]
        Angle [Degree] of the polarization filter in front of the detector.
    """
    coherent: bool = False
    mean_coupling: bool = False
    rotation: Quantity = field(default=0 * degree, init=False)
    mode_number: str = field(default='NC00', init=False)

    def __post_init__(self):
        self.binding = BindedDetector(
            mode_number=self.mode_number,
            sampling=self.sampling.magnitude,
            NA=self.NA.magnitude,
            cache_NA=self.cache_NA.magnitude,
            phi_offset=self.phi_offset.to(radian),
            gamma_offset=self.gamma_offset.to(radian),
            polarization_filter=self.polarization_filter.to(radian),
            rotation=self.rotation.to(radian),
            coherent=self.coherent,
            mean_coupling=self.mean_coupling
        )

    def get_structured_scalarfield(self, sampling: Optional[int] = 100) -> numpy.ndarray:
        """
        Generate a structured scalar field as a numpy array.

        Parameters
        ----------
            sampling : int
                The sampling rate for the scalar field. Default is 100.

        Returns
        -------
        numpy.ndarray
            A 2D array representing the structured scalar field.
        """
        return numpy.ones([sampling, sampling])


@dataclass(config=config_dict)
class IntegratingSphere(Photodiode):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This implies independence from the phase of the impinging scattered light field.

    Parameters
    ----------
    sampling: int)
        Sampling rate of the far-field distribution. Default is 200.
    polarization_filter: Union[float, None])
        Angle [Degree] of the polarization filter in front of the detector.
    """
    NA: Quantity = field(default=2 * AU, init=False)
    gamma_offset: Quantity = field(default=0 * degree, init=False)
    phi_offset: Quantity = field(default=0 * degree, init=False)

    def __post_init__(self):
        super(IntegratingSphere, self).__post_init__()

    def get_structured_scalarfield(self, sampling: Optional[int] = 100) -> numpy.ndarray:
        """
        Generate a structured scalar field as a numpy array.

        Parameters
        ----------
        sampling : int
            The sampling rate for the scalar field. Default is 100.

        Returns
        -------
        numpy.ndarray
            A 2D array representing the structured scalar field.
        """
        return numpy.ones([sampling, sampling])
