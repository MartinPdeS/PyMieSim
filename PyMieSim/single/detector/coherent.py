#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from typing import Optional
from pydantic.dataclasses import dataclass
from TypedUnit import Dimensionless, RefractiveIndex, Angle, ureg

from PyMieSim.utils import config_dict
from PyMieSim.binary.interface_detector import DETECTOR
from PyMieSim.single.detector.base import BaseDetector


@dataclass(config=config_dict, kw_only=True)
class CoherentMode(DETECTOR, BaseDetector):
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

    mode_number: str
    NA: Dimensionless
    gamma_offset: Angle
    phi_offset: Angle
    sampling: int = 200
    polarization_filter: Optional[Angle] = numpy.nan * ureg.degree
    cache_NA: Optional[Dimensionless] = 0 * ureg.AU
    mean_coupling: bool = False
    rotation: Angle = 90 * ureg.degree
    medium_refractive_index: RefractiveIndex = 1.0 * ureg.RIU

    def __post_init__(self):
        """
        Initialize the CoherentMode detector with its parameters.
        """

        if self.NA > 0.3 * ureg.AU or self.NA < 0 * ureg.AU:
            logging.warning(
                f"High values of NA: {self.NA} do not comply with the paraxial approximation. Values under 0.3 are preferred."
            )

        self.mode_family = self.mode_number[:2]

        if self.mode_family.lower() not in ["lp", "lg", "hg"]:
            raise ValueError(
                f'Invalid mode family: {self.mode_family}. Options are ["LP", "LG", "HG"]'
            )

        super().__init__(
            mode_number=self.mode_number,
            sampling=self.sampling,
            NA=self.NA.to(ureg.AU).magnitude,
            cache_NA=self.cache_NA.to(ureg.AU).magnitude,
            phi_offset=self.phi_offset.to(ureg.radian).magnitude,
            gamma_offset=self.gamma_offset.to(ureg.radian).magnitude,
            polarization_filter=self.polarization_filter.to(ureg.radian).magnitude,
            rotation=self.rotation.to(ureg.radian).magnitude,
            is_coherent=True,
            mean_coupling=self.mean_coupling,
            medium_refractive_index=self.medium_refractive_index.to(ureg.RIU).magnitude,
        )


# -
