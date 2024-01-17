#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from dataclasses import dataclass

from DataVisual import Xparameter

from PyMieSim import load_lp_mode
import PyMieSim.datavisual_x_parameters as Kwargs

from PyMieSim.binary.Sets import CppDetectorSet


@dataclass
class Photodiode():
    NA: list
    """ Numerical aperture of imaging system. """
    gamma_offset: list
    """ Angle [degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: list
    """ Angle [degree] offset of detector in the direction parallel to polarization. """
    polarization_filter: list
    """ Angle [degree] of polarization filter in front of detector. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = False
    """ Describe the detection scheme coherent or uncoherent. """
    sampling: int = 200
    """ Sampling number for the field evaluation. """
    name = "Photodiode"
    """ name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        self.initialize_binding()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.scalarfield = numpy.ones([1, self.sampling])

        self.NA = Xparameter(values=self.NA, **Kwargs.NA)

        self.phi_offset = Xparameter(values=self.phi_offset, **Kwargs.phi_offset)

        self.gamma_offset = Xparameter(values=self.gamma_offset, **Kwargs.gamma_offset)

        self.polarization_filter = Xparameter(values=self.polarization_filter, **Kwargs.polarization_filter)

        self.scalarfield = Xparameter(
            values=self.scalarfield,
            representation=numpy.asarray(['Photodiode']),
            **Kwargs.scalarfield
        )

    def initialize_binding(self):
        point_coupling = True if self.coupling_mode == 'point' else False
        phi_offset_rad = numpy.deg2rad(self.phi_offset.values)
        gamma_offset_rad = numpy.deg2rad(self.gamma_offset.values)

        self.binding = CppDetectorSet(
            scalarfield=self.scalarfield.values,
            NA=self.NA.values,
            phi_offset=phi_offset_rad,
            gamma_offset=gamma_offset_rad,
            polarization_filter=self.polarization_filter.values,
            point_coupling=point_coupling,
            coherent=self.coherent
        )

    def bind_to_experiment(self, experiment):
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_detector(self.binding)

    def append_to_table(self, table):
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, self.scalarfield, self.NA, self.phi_offset, self.gamma_offset, self.polarization_filter]

    def format_inputs(self):
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


@dataclass
class LPMode():
    mode_number: list
    """ List of mode to be used. """
    NA: list
    """ Numerical aperture of imaging system. """
    gamma_offset: list
    """ Angle [degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: list
    """ Angle [degree] offset of detector in the direction parallel to polarization. """
    polarization_filter: list
    """ Angle [degree] of polarization filter in front of detector. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = True
    """ Describe the detection scheme coherent or uncoherent. """
    sampling: int = 200
    """ Sampling number for the field evaluation. """
    name = "LPMode"
    """ name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        self.initialize_binding()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        representation = [
            f"{mode}" for mode in self.mode_number
        ]

        self.NA = Xparameter(values=self.NA, **Kwargs.NA)

        self.phi_offset = Xparameter(values=self.phi_offset, **Kwargs.phi_offset)

        self.gamma_offset = Xparameter(values=self.gamma_offset, **Kwargs.gamma_offset)

        self.polarization_filter = Xparameter(values=self.polarization_filter, **Kwargs.polarization_filter)

        self.scalarfield = Xparameter(
            values=self.get_fields_array(),
            representation=representation,
            **Kwargs.scalarfield
        )

    def get_fields_array(self) -> numpy.ndarray:
        """
        Gets the field arrays for the detection schemes the first dimension is the different fields.
        Second dimension are the individual mode fields.

        :returns:   The fields array.
        :rtype:     numpy.ndarray
        """
        fields_array = [
            load_lp_mode(
                mode_number=mode,
                sampling=self.sampling,
                structure_type='unstructured'
            ) for mode in self.mode_number
        ]

        return numpy.asarray(fields_array).astype(numpy.complex128)

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
            scalarfield=self.scalarfield.values,
            NA=self.NA.values,
            phi_offset=phi_offset_rad,
            gamma_offset=gamma_offset_rad,
            polarization_filter=self.polarization_filter.values,
            point_coupling=point_coupling,
            coherent=self.coherent
        )

    def bind_to_experiment(self, experiment) -> None:
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_detector(self.binding)

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, self.scalarfield, self.NA, self.phi_offset, self.gamma_offset, self.polarization_filter]

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)

        self.NA = numpy.atleast_1d(self.NA).astype(float)

        self.phi_offset = numpy.atleast_1d(self.phi_offset).astype(float)

        self.gamma_offset = numpy.atleast_1d(self.gamma_offset).astype(float)

        self.polarization_filter = numpy.atleast_1d(self.polarization_filter).astype(float)
