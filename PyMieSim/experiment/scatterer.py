#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.experiment.setup import Setup
    from PyMieSim.experiment.source import Gaussian, PlaneWave


from collections.abc import Iterable

import numpy
from dataclasses import dataclass

from DataVisual import Xparameter

import PyMieSim.datavisual_x_parameters as Kwargs

from PyMieSim.binary.Sets import CppCoreShellSet, CppCylinderSet, CppSphereSet


class BaseScatterer():
    def __post_init__(self):
        self.asserts_inputs()

        self.format_inputs()

        self.build_x_parameters()

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        for parameter_str in self.parameter_str_list:
            parameter = getattr(self, parameter_str)

            parameter = numpy.atleast_1d(parameter)

            setattr(self, parameter_str, parameter)

    def bind_to_experiment(self, experiment: Setup) -> None:
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        method_str = 'set_' + self.name

        getattr(experiment.binding, method_str)(self.binding)

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No returns
        :rtype:     None
        """
        for parameter_str in self.parameter_str_list:
            parameter = getattr(self, parameter_str)

            kwargs_parameter = getattr(Kwargs, parameter_str)

            x_parameter = Xparameter(values=parameter, **kwargs_parameter)

            setattr(self, parameter_str, x_parameter)


@dataclass
class CoreShell(BaseScatterer):
    core_diameter: Iterable
    """ diameter of the core of the single scatterer [m]. """
    shell_width: Iterable
    """ Width of the shell of the single scatterer [m]. """
    core_index: Iterable = None
    """ Refractive index of the core of the scatterer. """
    shell_index: Iterable = None
    """ Refractive index of the shell of the scatterer. """
    core_material: Iterable = None
    """ Core material of which the scatterer is made of. Only if core_index is not specified.  """
    shell_material: Iterable = None
    """ Shell material of which the scatterer is made of. Only if shell_index is not specified.  """
    n_medium: list = 1.0
    """ Refractive index of scatterer medium. """
    name: str = 'coreshell'
    """name of the set """

    parameter_str_list = [
        'n_medium',
        'core_diameter',
        'shell_width',
        'core_material',
        'shell_material',
        'core_index',
        'shell_index'
    ]

    def asserts_inputs(self) -> None:
        """
        Asserts that core and shell are either defined by index or material

        :returns:   No returns
        :rtype:     None
        """
        assert (self.core_material is None) ^ (self.core_index is None), "core_material xor core_index has to be defined"

        assert (self.shell_material is None) ^ (self.shell_index is None), "shell_material xor shell_index has to be defined"

        self.bounded_core = True if self.core_material is not None else False

        self.bounded_shell = True if self.shell_material is not None else False

    def get_core_shell_diameter_from_shell_width(
            self,
            core_diameter: numpy.ndarray,
            shell_width: numpy.ndarray) -> tuple:
        """
        Gets the core shell diameter from shell width.

        :param      core_diameter:  The core diameter
        :type       core_diameter:  numpy.ndarray
        :param      shell_width:    The shell width
        :type       shell_width:    numpy.ndarray

        :returns:   The core shell diameter from shell width.
        :rtype:     tuple
        """
        _core_diameter = numpy.tile(core_diameter, reps=len(shell_width))
        shell_diameter = numpy.add.outer(core_diameter, shell_width)
        return _core_diameter.flatten(), shell_diameter.T.flatten()

    def bind_core_material_shell_material(self, source: Gaussian | PlaneWave) -> None:
        core_material = [
            material.GetRI(source.values) for material in self.core_material
        ]

        shell_material = [
            material.GetRI(source.values) for material in self.shell_material
        ]

        core_material = numpy.asarray(core_material)
        shell_material = numpy.asarray(shell_material)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values.astype(float),
            shell_width=self.shell_width.values.astype(float),
            core_material=core_material.astype(numpy.complex128),
            shell_material=shell_material.astype(numpy.complex128),
            n_medium=self.n_medium.values.astype(float),
        )

    def bind_core_material_shell_index(self, source) -> None:
        core_material = [
            material.GetRI(source.values) for material in self.core_material
        ]

        core_material = numpy.asarray(core_material)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values.astype(float),
            shell_width=self.shell_width.values.astype(float),
            core_material=core_material.astype(numpy.complex128),
            shell_index=self.shell_index.values.astype(numpy.complex128),
            n_medium=self.n_medium.values.astype(float),
        )

    def bind_index_material_shell_material(self, source: Gaussian | PlaneWave) -> None:
        shell_material = [
            material.GetRI(source.values) for material in self.shell_material
        ]

        shell_material = numpy.asarray(shell_material)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values.astype(float),
            shell_width=self.shell_width.values.astype(float),
            core_index=self.core_index.values.astype(numpy.complex128),
            shell_material=shell_material.astype(numpy.complex128),
            n_medium=self.n_medium.values.astype(float),
        )

    def bind_index_material_shell_index(self, source: Gaussian | PlaneWave) -> None:
        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values.astype(float),
            shell_width=self.shell_width.values.astype(float),
            core_index=self.core_index.values.astype(numpy.complex128),
            shell_index=self.shell_index.values.astype(numpy.complex128),
            n_medium=self.n_medium.values.astype(float),
        )

    def evaluate_index_material(self, source: Gaussian | PlaneWave) -> None:
        if self.bounded_core and self.bounded_shell:
            self.bind_core_material_shell_material(source)

        if self.bounded_core and not self.bounded_shell:
            self.bind_core_material_shell_index(source)

        if not self.bounded_core and self.bounded_shell:
            self.bind_index_material_shell_material(source)

        if not self.bounded_core and not self.bounded_shell:
            self.bind_index_material_shell_index(source)

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        if self.bounded_core and self.bounded_shell:
            index_material_table = [self.core_material, self.shell_material]

        if self.bounded_core and not self.bounded_shell:
            index_material_table = [self.core_material, self.shell_index]

        if not self.bounded_core and self.bounded_shell:
            index_material_table = [self.core_index, self.shell_material]

        if not self.bounded_core and not self.bounded_shell:
            index_material_table = [self.core_index, self.shell_index]

        return [*table, self.core_diameter, self.shell_width, *index_material_table, self.n_medium]


@dataclass
class Sphere(BaseScatterer):
    diameter: Iterable
    """ diameter of the single scatterer in unit of meter. """
    index: Iterable = None
    """ Refractive index of scatterer. """
    material: Iterable = None
    """ material of which the scatterer is made of. Only if index is not specified. """
    n_medium: Iterable = 1.0
    """ Refractive index of scatterer medium. """
    name: str = 'sphere'
    """name of the set """

    parameter_str_list = [
        'n_medium',
        'diameter',
        'material',
        'index'
    ]

    def asserts_inputs(self) -> None:
        """
        Asserts that core and shell are either defined by index or material

        :returns:   No returns
        :rtype:     None
        """
        assert (self.material is None) ^ (self.index is None), "material xor index has to be defined"

        self.bounded_index = True if self.material is not None else False

    def evaluate_index_material(self, source: Gaussian | PlaneWave) -> None:
        if self.bounded_index:
            material = [
                material.GetRI(source.values) for material in self.material
            ]

            material_index = numpy.asarray(material).astype(complex)

            self.binding = CppSphereSet(
                diameter=self.diameter.values.astype(float),
                material=material_index.astype(numpy.complex128),
                n_medium=self.n_medium.values.astype(float)
            )

        else:
            self.binding = CppSphereSet(
                diameter=self.diameter.values.astype(float),
                index=self.index.values.astype(numpy.complex128),
                n_medium=self.n_medium.values.astype(float)
            )

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        if self.bounded_index:
            index_material_table = [self.material]

        else:
            index_material_table = [self.index]

        return [*table, self.diameter, *index_material_table, self.n_medium]


@dataclass
class Cylinder(BaseScatterer):
    diameter: Iterable
    """ diameter of the single scatterer in unit of meter. """
    index: Iterable = None
    """ Refractive index of scatterer. """
    material: Iterable = None
    """ Refractive index of scatterer medium. """
    n_medium: list = 1.0
    """ material of which the scatterer is made of. Only if index is not specified. """
    name: str = 'cylinder'
    """name of the set """

    parameter_str_list = [
        'n_medium',
        'diameter',
        'material',
        'index'
    ]

    def asserts_inputs(self) -> None:
        """
        Asserts that core and shell are either defined by index or material

        :returns:   No returns
        :rtype:     None
        """
        assert (self.material is None) ^ (self.index is None), "material xor index has to be defined"

        self.bounded_index = True if self.material is not None else False

    def evaluate_index_material(self, source: Gaussian | PlaneWave):
        if self.bounded_index:
            material = numpy.asarray(
                [
                    material.GetRI(source.values) for material in self.material
                ]
            )

            self.binding = CppCylinderSet(
                diameter=self.diameter.values.astype(float),
                material=material.astype(numpy.complex128),
                n_medium=self.n_medium.values.astype(float)
            )

        else:
            self.binding = CppCylinderSet(
                diameter=self.diameter.values.astype(float),
                index=self.index.values.astype(numpy.complex128),
                n_medium=self.n_medium.values.astype(float)
            )

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        if self.bounded_index:
            index_material_table = [self.material]

        else:
            index_material_table = [self.index]

        return [*table, self.diameter, *index_material_table, self.n_medium]

# -
