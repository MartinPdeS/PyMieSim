#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Optional
from PyMieSim.binary.interface_detector import DETECTOR
from PyMieSim.units import Quantity, degree, AU, radian, RIU
from PyMieSim.single.detector.base import BaseDetector


class Photodiode(DETECTOR, BaseDetector):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This means it is independent of the phase of the impinging scattered light field.

    """
    def __init__(self,
            NA: Quantity,
            gamma_offset: Quantity,
            phi_offset: Quantity,
            sampling: Quantity = 200,
            polarization_filter: Optional[Quantity] = None,
            cache_NA: Optional[Quantity] = 0 * AU,
            mean_coupling: bool = False,
            medium_refractive_index: Quantity = 1.0 * RIU,
        ):
        """
        Initialize the Photodiode detector with its parameters.

        Parameters
        ----------

        NA : Quantity
            Numerical aperture of the imaging system.
        gamma_offset : Quantity
            Angle [Degree] offset of the detector in the direction perpendicular to polarization.
        phi_offset : Quantity
            Angle [Degree] offset of the detector in the direction parallel to polarization.
        sampling : Quantity
            Sampling rate of the far-field distribution. Default is 200.
        polarization_filter : Optional[Quantity]
            Angle [Degree] of the polarization filter in front of the detector.
        cache_NA : Optional[Quantity]
            Numerical aperture of the detector cache. Default is 0 AU.
        mean_coupling : bool
            Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
        medium_refractive_index : Quantity
            The refractive index of the medium in which the detector operates. This is important for
            determining the acceptance cone of light.
            Default is 1.0 (vacuum or air).
        """

        self.NA = self._validate_AU_units(NA)
        self.cache_NA = self._validate_AU_units(cache_NA)
        self.gamma_offset = self._validate_angle_units(gamma_offset)
        self.phi_offset = self._validate_angle_units(phi_offset)
        self.sampling = self._validate_no_units(sampling)
        self.polarization_filter = self._validate_polarization_units(polarization_filter)
        self.mean_coupling = mean_coupling
        self.medium_refractive_index = self._validate_RIU_units(medium_refractive_index)

        super().__init__(
            mode_number="NC00",
            NA=self.NA.to_base_units().magnitude,
            cache_NA=self.cache_NA.to_base_units().magnitude,
            gamma_offset=self.gamma_offset.to(radian).magnitude,
            phi_offset=self.phi_offset.to(radian).magnitude,
            sampling=self.sampling,
            polarization_filter=self.polarization_filter.to(radian).magnitude,
            mean_coupling=self.mean_coupling,
            rotation=0,
            coherent=False,
            medium_refractive_index=self.medium_refractive_index.to_base_units().magnitude
        )


class IntegratingSphere(Photodiode):
    def __init__(self, sampling: Quantity, polarization_filter: Optional[Quantity] = None, mean_coupling: bool = False):
        """
        Detector class representing a photodiode with a non-coherent light coupling mechanism.
        This implies independence from the phase of the impinging scattered light field.

        Parameters
        ----------

        sampling : Quantity
            Sampling rate of the far-field distribution. Default is 200.
        polarization_filter : Optional[Quantity]
            Angle [Degree] of the polarization filter in front of the detector.
        cache_NA : Optional[Quantity]
            Numerical aperture of the detector cache. Default is 0 AU.
        mean_coupling : bool
            Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
        """

        super().__init__(
            sampling=sampling,
            polarization_filter=polarization_filter,
            mean_coupling=mean_coupling,
            NA=2.0 * AU,  # Fixed NA for IntegratingSphere
            gamma_offset=0 * degree,  # Fixed gamma offset for IntegratingSphere
            phi_offset=0 * degree,  # Fixed phi offset for IntegratingSphere
            cache_NA=0 * AU,  # Fixed cache NA for IntegratingSphere
        )
