#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from typing import Optional
from TypedUnit import RefractiveIndex, Dimensionless, Angle, ureg

from PyMieSim.binary.interface_single import DETECTOR


class Photodiode(DETECTOR):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This implies independence from the phase of the impinging scattered light field.
    """
    def __init__(
        self, NA: Dimensionless,
        gamma_offset: Angle,
        phi_offset: Angle,
        medium_refractive_index: RefractiveIndex,
        sampling: int = 200,
        polarization_filter: Optional[Angle] = numpy.nan * ureg.degree,
        cache_NA: Optional[Dimensionless] = 0 * ureg.AU,
        mean_coupling: bool = False,
    ):
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
        """
        NA = Dimensionless.check(NA)
        gamma_offset = Angle.check(gamma_offset)
        phi_offset = Angle.check(phi_offset)
        sampling = int(sampling)

        if polarization_filter is not None:
            polarization_filter = Angle.check(polarization_filter)
        else:
            polarization_filter = numpy.nan * ureg.degree

        if cache_NA is not None:
            cache_NA = Dimensionless.check(cache_NA)
        else:
            cache_NA = 0 * ureg.AU

        medium_refractive_index = RefractiveIndex.check(medium_refractive_index)

        super().__init__(
            mode_number="NC00",
            NA=NA,
            cache_NA=cache_NA,
            gamma_offset=gamma_offset,
            phi_offset=phi_offset,
            sampling=sampling,
            polarization_filter=polarization_filter,
            mean_coupling=mean_coupling,
            rotation=0 * ureg.degree,
            is_coherent=False,
            medium_refractive_index=medium_refractive_index,
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
            medium_refractive_index=1.8 * ureg.RIU
        )
