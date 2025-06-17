#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Optional
from PyMieSim.binary.interface_detector import DETECTOR
from PyMieSim import units
from PyMieSim.single.detector.base import BaseDetector


class Photodiode(DETECTOR, BaseDetector):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This means it is independent of the phase of the impinging scattered light field.

    """
    def __init__(self,
            NA: units.Quantity,
            gamma_offset: units.Quantity,
            phi_offset: units.Quantity,
            sampling: units.Quantity = 200,
            polarization_filter: Optional[units.Quantity] = None,
            cache_NA: Optional[units.Quantity] = 0 * units.AU,
            mean_coupling: bool = False,
            medium_refractive_index: units.Quantity = 1.0 * units.RIU):
        """
        Initialize the Photodiode detector with its parameters.

        Parameters
        ----------

        NA : units.Quantity
            Numerical aperture of the imaging system.
        gamma_offset : units.Quantity
            Angle [Degree] offset of the detector in the direction perpendicular to polarization.
        phi_offset : units.Quantity
            Angle [Degree] offset of the detector in the direction parallel to polarization.
        sampling : units.Quantity
            Sampling rate of the far-field distribution. Default is 200.
        polarization_filter : Optional[units.Quantity]
            Angle [Degree] of the polarization filter in front of the detector.
        cache_NA : Optional[units.Quantity]
            Numerical aperture of the detector cache. Default is 0 AU.
        mean_coupling : bool
            Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
        medium_refractive_index : units.Quantity
            The refractive index of the medium in which the detector operates. This is important for
            determining the acceptance cone of light.
            Default is 1.0 (vacuum or air).
        """
        self.mean_coupling = mean_coupling

        self.NA = self._validate_units(NA, dimension='arbitrary', units=units.AU)
        self.cache_NA = self._validate_units(cache_NA, dimension='arbitrary', units=units.AU)
        self.sampling = self._validate_units(sampling, dimension='arbitrary', units=units.AU)

        self.gamma_offset = self._validate_units(gamma_offset, dimension='angle', units=units.degree)
        self.phi_offset = self._validate_units(phi_offset, dimension='angle', units=units.degree)

        self.medium_refractive_index = self._validate_units(medium_refractive_index, dimension='refractive index', units=units.RIU)

        self.polarization_filter = self._validate_detector_polarization_units(polarization_filter)

        super().__init__(
            mode_number="NC00",
            NA=self.NA.to(units.AU).magnitude,
            cache_NA=self.cache_NA.to(units.AU).magnitude,
            gamma_offset=self.gamma_offset.to(units.radian).magnitude,
            phi_offset=self.phi_offset.to(units.radian).magnitude,
            sampling=self.sampling,
            polarization_filter=self.polarization_filter.to(units.radian).magnitude,
            mean_coupling=self.mean_coupling,
            rotation=0,
            coherent=False,
            medium_refractive_index=self.medium_refractive_index.to(units.RIU).magnitude
        )


class IntegratingSphere(Photodiode):
    def __init__(self, sampling: units.Quantity, polarization_filter: Optional[units.Quantity] = None, mean_coupling: bool = False):
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
            NA=2.0 * units.AU,  # Fixed NA for IntegratingSphere
            gamma_offset=0 * units.degree,  # Fixed gamma offset for IntegratingSphere
            phi_offset=0 * units.degree,  # Fixed phi offset for IntegratingSphere
            cache_NA=0 * units.AU,  # Fixed cache NA for IntegratingSphere
        )
