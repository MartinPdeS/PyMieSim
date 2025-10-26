#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from typing import Optional
from pydantic.dataclasses import dataclass
from TypedUnit import RefractiveIndex, Dimensionless, Angle, ureg

from PyMieSim.utils import config_dict
from PyMieSim.binary.interface_detector import DETECTOR
from PyMieSim.single.detector.base import BaseDetector


@dataclass(config=config_dict, kw_only=True)
class Photodiode(DETECTOR, BaseDetector):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This means it is independent of the phase of the impinging scattered light field.

    Parameters
    ----------

    NA : Dimensionless
        Numerical aperture of the imaging system.
    gamma_offset : Angle
        Angle [Degree] offset of the detector in the direction perpendicular to polarization.
    phi_offset : Angle
        Angle [Degree] offset of the detector in the direction parallel to polarization.
    sampling : int
        Sampling rate of the far-field distribution. Default is 200.
    polarization_filter : Optional[Angle]
        Angle [Degree] of the polarization filter in front of the detector.
    cache_NA : Optional[Dimensionless]
        Numerical aperture of the detector cache. Default is 0 AU.
    mean_coupling : bool
        Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
    medium_refractive_index : RefractiveIndex
        The refractive index of the medium in which the detector operates. This is important for
        determining the acceptance cone of light.
        Default is 1.0 (vacuum or air).
    """

    NA: Dimensionless
    gamma_offset: Angle
    phi_offset: Angle
    sampling: int = 200
    polarization_filter: Optional[Angle] = numpy.nan * ureg.degree
    cache_NA: Optional[Dimensionless] = 0 * ureg.AU
    mean_coupling: bool = False
    medium_refractive_index: RefractiveIndex = 1.0 * ureg.RIU

    def __post_init__(self):
        """
        Initialize the Photodiode detector with its parameters.
        """
        super().__init__(
            mode_number="NC00",
            NA=self.NA.to(ureg.AU).magnitude,
            cache_NA=self.cache_NA.to(ureg.AU).magnitude,
            gamma_offset=self.gamma_offset.to(ureg.radian).magnitude,
            phi_offset=self.phi_offset.to(ureg.radian).magnitude,
            sampling=self.sampling,
            polarization_filter=self.polarization_filter.to(ureg.radian).magnitude,
            mean_coupling=self.mean_coupling,
            rotation=0,
            is_coherent=False,
            medium_refractive_index=self.medium_refractive_index.to(ureg.RIU).magnitude,
        )


class IntegratingSphere(Photodiode):
    def __init__(
        self,
        sampling: ureg.Quantity,
        polarization_filter: Optional[Angle] = numpy.nan * ureg.degree,
        mean_coupling: bool = False,
    ):
        """
        Detector class representing a photodiode with a non-coherent light coupling mechanism.
        This implies independence from the phase of the impinging scattered light field.

        Parameters
        ----------

        sampling : units.Quantity
            Sampling rate of the far-field distribution. Default is 200.
        polarization_filter : Optional[units.Quantity]
            Angle [Degree] of the polarization filter in front of the detector.
        cache_NA : Optional[units.Quantity]
            Numerical aperture of the detector cache. Default is 0 AU.
        mean_coupling : bool
            Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
        """

        super().__init__(
            sampling=sampling,
            polarization_filter=polarization_filter,
            mean_coupling=mean_coupling,
            NA=2.0 * ureg.AU,  # Fixed NA for IntegratingSphere
            gamma_offset=0 * ureg.degree,  # Fixed gamma offset for IntegratingSphere
            phi_offset=0 * ureg.degree,  # Fixed phi offset for IntegratingSphere
            cache_NA=0 * ureg.AU,  # Fixed cache NA for IntegratingSphere
        )
