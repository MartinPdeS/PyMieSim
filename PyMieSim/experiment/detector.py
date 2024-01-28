#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.experiment.setup import Setup
    from collections.abc import Iterable

import numpy
from dataclasses import dataclass, field

from DataVisual import Xparameter
from PyMieSim import load_lp_mode
import PyMieSim.datavisual_x_parameters as Kwargs
from PyMieSim.binary.Sets import CppDetectorSet


@dataclass
class BaseDetector():
    NA: Iterable
    """ Numerical aperture of imaging system. """
    gamma_offset: Iterable
    """ Angle [degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: Iterable
    """ Angle [degree] offset of detector in the direction parallel to polarization. """
    polarization_filter: Iterable
    """ Angle [degree] of polarization filter in front of detector. """
    sampling: int
    """ Sampling number for the field evaluation. """

    parameter_str_list = [
        'NA',
        'phi_offset',
        'gamma_offset',
        'polarization_filter',
    ]

    def __post_init__(self):
        self.format_inputs()

        self.get_fields_array()

        self.get_rotation_angle_from_mode_number()

        self.build_x_parameters()

        self.initialize_binding()

    def get_rotation_angle_from_mode_number(self) -> None:
        """
        Compute the rotation angle from mode number if detector is LPMode else
        rotation angle is zero.

        :returns:   No return
        :rtype:     None
        """
        if self.name.lower() == 'photodiode':
            self.rotation_angle = numpy.zeros(self.scalarfield.values.size)
            return

        rotation_angle_list = []

        self.mode_number = [
            mode if ":" in mode else mode + ":00" for mode in self.mode_number
        ]

        for mode_number in self.mode_number:
            _, rotation_angle = mode_number.split(':')

            rotation_angle_list.append(rotation_angle)

        self.rotation_angle = numpy.asarray(rotation_angle_list).astype(float)

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.NA = numpy.atleast_1d(self.NA).astype(float)
        self.phi_offset = numpy.atleast_1d(self.phi_offset).astype(float)
        self.gamma_offset = numpy.atleast_1d(self.gamma_offset).astype(float)
        self.polarization_filter = numpy.atleast_1d(self.polarization_filter).astype(float)

    def bind_to_experiment(self, experiment: Setup) -> None:
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_detector(self.binding)

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No return
        :rtype:     None
        """
        for parameter_str in self.parameter_str_list:
            parameter = getattr(self, parameter_str)

            parameter = numpy.array(parameter)

            kwargs_parameter = getattr(Kwargs, parameter_str)

            x_parameter = Xparameter(values=parameter, **kwargs_parameter)

            setattr(self, parameter_str, x_parameter)

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, self.scalarfield, self.NA, self.phi_offset, self.gamma_offset, self.polarization_filter]

    def initialize_binding(self) -> None:
        """
        Initializes the cpp binding of the LPmode detector set.

        :returns:   No return
        :rtype:     None
        """
        point_coupling = True if self.coupling_mode == 'point' else False

        phi_offset_rad = numpy.deg2rad(self.phi_offset.values)

        gamma_offset_rad = numpy.deg2rad(self.gamma_offset.values)

        self.binding = CppDetectorSet(
            scalarfield=self.scalarfield.values.astype(complex),
            NA=self.NA.values,
            phi_offset=phi_offset_rad,
            gamma_offset=gamma_offset_rad,
            polarization_filter=self.polarization_filter.values,
            point_coupling=point_coupling,
            coherent=self.coherent,
            rotation_angle=self.rotation_angle.astype(float)
        )


@dataclass
class Photodiode(BaseDetector):
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = field(default=False, init=False)
    """ Describe the detection scheme coherent or uncoherent. """
    name: str = field(default="Photodiode", init=False)
    """ name of the set """

    def get_fields_array(self) -> numpy.ndarray:
        """
        Gets the field arrays for the detection schemes the first dimension is the different fields.
        Second dimension are the individual mode fields.

        :returns:   The fields array.
        :rtype:     numpy.ndarray
        """
        scalarfield = numpy.ones([1, self.sampling])

        self.scalarfield = Xparameter(
            values=scalarfield,
            representation=['Photodiode'],
            **Kwargs.scalarfield
        )


@dataclass
class LPMode(BaseDetector):
    mode_number: Iterable
    """ List of mode to be used. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = field(default=True, init=False)
    """ Describe the detection scheme coherent or uncoherent. """
    name: str = field(default="LPMode", init=False)
    """ name of the set """

    def get_fields_array(self) -> numpy.ndarray:
        """
        Gets the field arrays for the detection schemes the first dimension is the different fields.
        Second dimension are the individual mode fields.

        :returns:   The fields array.
        :rtype:     numpy.ndarray
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)

        scalarfield = load_lp_mode(
            mode_number=self.mode_number,
            sampling=self.sampling,
            structure_type='unstructured'
        ).astype(complex)

        self.scalarfield = Xparameter(
            values=scalarfield,
            representation=self.mode_number,
            **Kwargs.scalarfield
        )
