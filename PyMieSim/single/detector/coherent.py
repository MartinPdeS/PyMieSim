#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from typing import Optional
from pydantic.dataclasses import dataclass
from PyMieSim.binary.interface_detector import BindedDetector # type: ignore
from PyMieSim.binary import interface_mode_field
from PyMieSim.units import radian
from PyMieSim.single.detector.base import BaseDetector, config_dict


@dataclass(config=config_dict)
class CoherentMode(BaseDetector):
    """
    Detector class representing a laser Hermite-Gauss mode with a coherent light coupling mechanism.
    This means it depends on the phase of the impinging scattered light field.

    Parameters
    ----------
    mode_number: str
        String representing the HG mode to be initialized (e.g., 'LP01', 'HG11', 'LG22').
    NA: float
        Numerical aperture of the imaging system.
    cache_NA : Quantity
        Numerical aperture of the detector cache.
    gamma_offset: float
        Angle [Degree] offset of the detector in the direction perpendicular to polarization.
    phi_offset: float
        Angle [Degree] offset of the detector in the direction parallel to polarization.
    sampling: int
        Sampling rate of the far-field distribution. Default is 200.
    polarization_filter: Union[float, None]
        Angle [Degree] of the polarization filter in front of the detector.
    mean_coupling: bool
        Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
    coherent: bool
        Indicates if the coupling mechanism is coherent. Default is True.
    rotation: float
        Rotation angle of the field along the axis of propagation. Default is 90.
    """
    mode_number: str
    coherent: bool = True
    mean_coupling: bool = True

    def __post_init__(self):
        if self.NA > 0.3 or self.NA < 0:
            logging.warning(f"High values of NA: {self.NA} do not comply with the paraxial approximation. Values under 0.3 are preferred.")

        self.mode_family = self.mode_number[:2]

        if self.mode_family.lower() not in ['lp', 'lg', 'hg']:
            raise ValueError(f'Invalid mode family: {self.mode_family}. Options are ["LP", "LG", "HG"]')

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

        self.binding = BindedDetector(
            mode_number=self.mode_number,
            sampling=self.sampling.magnitude,
            NA=self.NA.magnitude,
            cache_NA=self.cache_NA.magnitude,
            phi_offset=self.phi_offset.to(radian).magnitude,
            gamma_offset=self.gamma_offset.to(radian).magnitude,
            polarization_filter=self.polarization_filter.to(radian).magnitude,
            rotation=self.rotation.to(radian).magnitude,
            coherent=True,
            mean_coupling=self.mean_coupling,
            medium_refractive_index=self.medium_refractive_index.magnitude
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
