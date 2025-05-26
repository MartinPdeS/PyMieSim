#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from typing import Optional
from PyMieSim.binary.interface_detector import DETECTOR
from PyMieSim.binary import interface_mode_field
from PyMieSim.units import radian, Quantity, degree, AU, RIU
from PyMieSim.single.detector.base import BaseDetector


class CoherentMode(DETECTOR, BaseDetector):
    def __init__(self,
            mode_number: str,
            NA: Quantity,
            gamma_offset: Quantity,
            phi_offset: Quantity,
            sampling: Quantity = 200 * AU,
            polarization_filter: Optional[Quantity] = None,
            cache_NA: Optional[Quantity] = 0 * AU,
            mean_coupling: bool = False,
            rotation: Quantity = 90 * degree,
            medium_refractive_index: Quantity = 1.0 * RIU,
        ):
        """
        Initialize the CoherentMode detector with its parameters.

        Depending on the mode_number, this detector can represent different types of modes such as Hermite-Gauss (HG), Laguerre-Gauss (LG), or fiber linearly polarized (LP).
        This class is designed to handle coherent light coupling mechanisms, meaning it depends on the phase of the impinging scattered light field.

        Parameters
        ----------
        mode_number : str
            String representing the mode to be initialized (e.g., 'LP01', 'HG11', 'LG22').
        NA : Quantity
            The numerical aperture of the detector, a dimensionless number that defines its light-gathering
            ability. Higher NA values indicate a wider acceptance angle for capturing light.
        gamma_offset : Quantity
            The rotational offset in the gamma direction, controlling the angular positioning of the detector relative to the scatterer.
        phi_offset : Quantity
            The rotational offset in the phi direction, specifying the azimuthal angle of the detector's orientation.
        sampling : Quantity
            The sampling resolution of the detector, controlling how finely the detector's field is sampled
            over the surface. A higher value increases precision but may require more computational resources. Default is 200.
        polarization_filter : Optional[Quantity]
            The polarization filter applied to the detected light. If specified, it allows the detector to
            selectively capture light with a particular polarization. If set to `nan`, no filtering is applied.
        cache_NA : Optional[Quantity]
            Numerical aperture of the detector cache. Default is 0 AU.
        mean_coupling : bool
            A flag indicating whether the mean value of the coupling is used when calculating the light
            interaction with the scatterer. This is typically set for cases where an average coupling is
            needed over multiple angles or configurations.
        rotation : Quantity
            The rotational angle of the detector, defining its orientation relative to the incoming light or
            scatterer. The default value rotates the detector by 90 degrees.
        medium_refractive_index : Quantity
            The refractive index of the medium in which the detector operates. This is important for
            determining the acceptance cone of light.
            Default is 1.0 (vacuum or air).
        """

        if NA > 0.3 * AU or NA < 0 * AU:
            logging.warning(f"High values of NA: {NA} do not comply with the paraxial approximation. Values under 0.3 are preferred.")

        self.mode_family = mode_number[:2]

        if self.mode_family.lower() not in ['lp', 'lg', 'hg']:
            raise ValueError(f'Invalid mode family: {self.mode_family}. Options are ["LP", "LG", "HG"]')

        self.mode_number = mode_number
        self.NA = self._validate_UA_units(NA)
        self.cache_NA = self._validate_UA_units(cache_NA)
        self.gamma_offset = self._validate_angle_units(gamma_offset)
        self.phi_offset = self._validate_angle_units(phi_offset)
        self.sampling = self._validate_UA_units(sampling)
        self.polarization_filter = self._validate_polarization_units(polarization_filter)
        self.mean_coupling = mean_coupling
        self.rotation = self._validate_angle_units(rotation)
        self.medium_refractive_index = self._validate_RIU_units(medium_refractive_index)

        super().__init__(
            mode_number=self.mode_number,
            NA=self.NA.to_base_units().magnitude,
            cache_NA=self.cache_NA.to_base_units().magnitude,
            gamma_offset=self.gamma_offset.to(radian).magnitude,
            phi_offset=self.phi_offset.to(radian).magnitude,
            sampling=self.sampling.to_base_units().magnitude,
            polarization_filter=self.polarization_filter.to(radian).magnitude,
            mean_coupling=self.mean_coupling,
            rotation=self.rotation.to(radian).magnitude,
            coherent=False,
            medium_refractive_index=self.medium_refractive_index.to_base_units().magnitude
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
        number_0, number_1 = self.mode_number[2:]
        self.number_0, self.number_1 = int(number_0), int(number_1)

        match self.mode_family.lower():
            case 'lp':
                self.azimuthal_number, self.radial_number = self.number_0, self.number_1
                self.cpp_mode_field_getter = interface_mode_field.get_LP
            case 'lg':
                self.azimuthal_number, self.radial_number = self.number_0, self.number_1
                self.cpp_mode_field_getter = interface_mode_field.get_LG
            case 'hg':
                self.x_number, self.y_number = self.number_0, self.number_1
                self.cpp_mode_field_getter = interface_mode_field.get_HG


        x_mesh, y_mesh = numpy.mgrid[-100:100:complex(sampling), -100:100:complex(sampling)]

        coordinates = numpy.row_stack((
            x_mesh.ravel(),
            y_mesh.ravel(),
        ))

        norm = numpy.sqrt(numpy.square(coordinates).sum(axis=0)).max()

        coordinates /= norm

        field = self.cpp_mode_field_getter(
            coordinates[0],
            coordinates[1],
            self.number_0,
            self.number_1
        )

        return field.reshape([sampling, sampling])


# -
