#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from typing import Optional
from PyMieSim import units
from PyMieSim.binary.interface_detector import DETECTOR
from PyMieSim.single.detector.base import BaseDetector


class CoherentMode(DETECTOR, BaseDetector):
    def __init__(self,
            mode_number: str,
            NA: units.Quantity,
            gamma_offset: units.Quantity,
            phi_offset: units.Quantity,
            sampling: units.Quantity = 200,
            polarization_filter: Optional[units.Quantity] = None,
            cache_NA: Optional[units.Quantity] = 0 * units.AU,
            mean_coupling: bool = False,
            rotation: units.Quantity = 90 * units.degree,
            medium_refractive_index: units.Quantity = 1.0 * units.RIU,
        ):
        """
        Initialize the CoherentMode detector with its parameters.

        Depending on the mode_number, this detector can represent different types of modes such as Hermite-Gauss (HG), Laguerre-Gauss (LG), or fiber linearly polarized (LP).
        This class is designed to handle coherent light coupling mechanisms, meaning it depends on the phase of the impinging scattered light field.

        Parameters
        ----------
        mode_number : str
            String representing the mode to be initialized (e.g., 'LP01', 'HG11', 'LG22').
        NA : units.Quantity
            The numerical aperture of the detector, a dimensionless number that defines its light-gathering
            ability. Higher NA values indicate a wider acceptance angle for capturing light.
        gamma_offset : units.Quantity
            The rotational offset in the gamma direction, controlling the angular positioning of the detector relative to the scatterer.
        phi_offset : units.Quantity
            The rotational offset in the phi direction, specifying the azimuthal angle of the detector's orientation.
        sampling : units.Quantity
            The sampling resolution of the detector, controlling how finely the detector's field is sampled
            over the surface. A higher value increases precision but may require more computational resources. Default is 200.
        polarization_filter : Optional[units.Quantity]
            The polarization filter applied to the detected light. If specified, it allows the detector to
            selectively capture light with a particular polarization. If set to `nan`, no filtering is applied.
        cache_NA : Optional[units.Quantity]
            Numerical aperture of the detector cache. Default is 0 AU.
        mean_coupling : bool
            A flag indicating whether the mean value of the coupling is used when calculating the light
            interaction with the scatterer. This is typically set for cases where an average coupling is
            needed over multiple angles or configurations.
        rotation : units.Quantity
            The rotational angle of the detector, defining its orientation relative to the incoming light or
            scatterer. The default value rotates the detector by 90 degrees.
        medium_refractive_index : units.Quantity
            The refractive index of the medium in which the detector operates. This is important for
            determining the acceptance cone of light.
            Default is 1.0 (vacuum or air).
        """

        if NA > 0.3 * units.AU or NA < 0 * units.AU:
            logging.warning(f"High values of NA: {NA} do not comply with the paraxial approximation. Values under 0.3 are preferred.")

        self.mode_family = mode_number[:2]

        if self.mode_family.lower() not in ['lp', 'lg', 'hg']:
            raise ValueError(f'Invalid mode family: {self.mode_family}. Options are ["LP", "LG", "HG"]')

        self.mode_number = mode_number
        self.mean_coupling = mean_coupling

        self.NA = self._validate_units(NA, dimension='arbitrary', units=units.AU)
        self.cache_NA = self._validate_units(cache_NA, dimension='arbitrary', units=units.AU)
        self.sampling = self._validate_units(sampling, dimension='arbitrary', units=units.AU)

        self.gamma_offset = self._validate_units(gamma_offset, dimension='angle', units=units.degree)
        self.phi_offset = self._validate_units(phi_offset, dimension='angle', units=units.degree)
        self.rotation = self._validate_units(rotation, dimension='angle', units=units.degree)

        self.medium_refractive_index = self._validate_units(medium_refractive_index, dimension='refractive index', units=units.RIU)

        self.polarization_filter = self._validate_detector_polarization_units(polarization_filter)

        super().__init__(
            mode_number=self.mode_number,
            sampling=self.sampling,
            NA=self.NA.to(units.AU).magnitude,
            cache_NA=self.cache_NA.to(units.AU).magnitude,
            phi_offset=self.phi_offset.to(units.radian).magnitude,
            gamma_offset=self.gamma_offset.to(units.radian).magnitude,
            polarization_filter=self.polarization_filter.to(units.radian).magnitude,
            rotation=self.rotation.to(units.radian).magnitude,
            coherent=False,
            mean_coupling=self.mean_coupling,
            medium_refractive_index=self.medium_refractive_index.to(units.RIU).magnitude,
        )

# -
