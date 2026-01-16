#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from typing import Optional
from PyMieSim.units import Dimensionless, RefractiveIndex, Angle, ureg

from PyMieSim.binary.interface_single import DETECTOR


class CoherentMode(DETECTOR):
    """
    Detector class representing a coherent mode detector.

    Depending on the mode_number, this detector can represent different types of modes such as Hermite-Gauss (HG), Laguerre-Gauss (LG), or fiber linearly polarized (LP).
    This class is designed to handle coherent light coupling mechanisms, meaning it depends on the phase of the impinging scattered light field.

    Parameters
    ----------
    mode_number : str
        String representing the mode to be initialized (e.g., 'LP01', 'HG11', 'LG22').
    NA : Dimensionless
        The numerical aperture of the detector, a dimensionless number that defines its light-gathering
        ability. Higher NA values indicate a wider acceptance angle for capturing light.
    gamma_offset : Angle
        The rotational offset in the gamma direction, controlling the angular positioning of the detector relative to the scatterer.
    phi_offset : Angle
        The rotational offset in the phi direction, specifying the azimuthal angle of the detector's orientation.
    sampling : int
        The sampling resolution of the detector, controlling how finely the detector's field is sampled
        over the surface. A higher value increases precision but may require more computational resources. Default is 200.
    polarization_filter : Optional[Angle]
        The polarization filter applied to the detected light. If specified, it allows the detector to
        selectively capture light with a particular polarization. If set to `nan`, no filtering is applied.
    cache_NA : Optional[Dimensionless]
        Numerical aperture of the detector cache. Default is 0 AU.
    mean_coupling : bool
        A flag indicating whether the mean value of the coupling is used when calculating the light
        interaction with the scatterer. This is typically set for cases where an average coupling is
        needed over multiple angles or configurations.
    rotation : Angle
        The rotational angle of the detector, defining its orientation relative to the incoming light or
        scatterer. The default value rotates the detector by 90 degrees.
    medium_refractive_index : RefractiveIndex
        The refractive index of the medium in which the detector operates. This is important for
        determining the acceptance cone of light.
        Default is 1.0 (vacuum or air).
    """

    def __init__(
        self,
        mode_number: str,
        NA: Dimensionless,
        gamma_offset: Angle,
        phi_offset: Angle,
        sampling: int = 200,
        polarization_filter: Optional[Angle] = numpy.nan * ureg.degree,
        cache_NA: Optional[Dimensionless] = 0 * ureg.AU,
        mean_coupling: bool = False,
        rotation: Angle = 90 * ureg.degree,
        medium_refractive_index: RefractiveIndex = 1.0 * ureg.RIU,
    ):
        """
        Initialize the CoherentMode detector with its parameters.
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

        if NA > 0.3 * ureg.AU or NA < 0 * ureg.AU:
            logging.warning(
                f"High values of NA: {NA} do not comply with the paraxial approximation. Values under 0.3 are preferred."
            )

        self.mode_family = mode_number[:2]

        if self.mode_family.lower() not in ["lp", "lg", "hg"]:
            raise ValueError(
                f'Invalid mode family: {self.mode_family}. Options are ["LP", "LG", "HG"]'
            )

        super().__init__(
            mode_number=mode_number,
            sampling=sampling,
            NA=NA,
            cache_NA=cache_NA,
            phi_offset=phi_offset,
            gamma_offset=gamma_offset,
            polarization_filter=polarization_filter,
            rotation=rotation,
            is_coherent=True,
            mean_coupling=mean_coupling,
            medium_refractive_index=medium_refractive_index,
        )
