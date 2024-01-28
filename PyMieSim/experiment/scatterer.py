#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.experiment.setup import Setup
    from PyMieSim.experiment.source import Gaussian, PlaneWave
    from collections.abc import Iterable

import numpy
from dataclasses import dataclass, field

from DataVisual import Xparameter
import PyMieSim.datavisual_x_parameters as Kwargs
from PyMieSim.binary.Sets import CppCoreShellSet, CppCylinderSet, CppSphereSet


@dataclass
class BaseScatterer():
    n_medium: Iterable
    """ material of which the scatterer is made of. Only if index is not specified. """
    source_set: Gaussian | PlaneWave
    """ Source set for the scatterer properties to be evaluated. """

    def __post_init__(self):
        self.asserts_inputs()

        self.format_inputs()

        self.build_binding_kwargs()

        self.build_x_parameters()

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        for parameter_str in self.parameter_dictionnary.keys():
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

    def build_binding_kwargs(self) -> None:
        self.evaluate_index_material()

        self.binding_kwargs = dict()

        for parameter_str in self.cpp_binding_str:
            values = getattr(self, parameter_str)

            self.binding_kwargs[parameter_str] = values

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No return
        :rtype:     None
        """
        self.x_table = []
        for parameter_str, dic in self.parameter_dictionnary.items():
            values = getattr(self, parameter_str)

            values = numpy.asarray(values)

            kwargs_parameter = getattr(Kwargs, parameter_str)

            x_parameter = Xparameter(values=values, **kwargs_parameter)

            setattr(self, parameter_str, x_parameter)

            self.x_table.append(x_parameter)

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, *self.x_table]

    def initialize_binding(self) -> None:
        """
        Initializes the cpp binding of the scatterer.

        :returns:   No return
        :rtype:     None
        """
        self.binding = self.binding_class(**self.binding_kwargs)

    def evaluate_index_material(self) -> None:
        """
        Evaluates all sets parameters that finishes with "material"
        to create the material_index paramter

        :returns:   No return
        :rtype:     None
        """
        for parameter_str, value in self.parameter_dictionnary.items():
            if parameter_str.endswith('material'):
                materials = getattr(self, parameter_str)

                material_index = [
                    material.GetRI(self.source_set.wavelength.values) for material in materials
                ]

                material_index = numpy.asarray(material_index).astype(complex)

                setattr(self, parameter_str + '_index', material_index)


@dataclass
class Sphere(BaseScatterer):
    diameter: Iterable
    """ diameter of the single scatterer in unit of meter. """
    index: Iterable = None
    """ Refractive index of scatterer. """
    material: Iterable = None
    """ material of which the scatterer is made of. Only if index is not specified. """
    name: str = field(default="sphere", init=False)
    """ Name of the set """

    def __post_init__(self):
        self.cpp_binding_str: list = [
            'diameter',
            'material_index',
            'index',
            'n_medium',
        ]

        self.parameter_dictionnary: dict = dict(
            diameter=dict(type=float, value=None),
            material=dict(type=complex, value=None),
            index=dict(type=complex, value=None),
            n_medium=dict(type=float, value=None),
        )

        self.binding_class: type = CppSphereSet

        super().__post_init__()

    def asserts_inputs(self) -> None:
        """
        Asserts that core and shell are either defined by index or material

        :returns:   No return
        :rtype:     None
        """
        assert (self.material is None) ^ (self.index is None), "material xor index has to be defined"

        if self.material is not None:
            del self.parameter_dictionnary['index']
            self.cpp_binding_str.remove('index')
        else:
            del self.parameter_dictionnary['material']
            self.cpp_binding_str.remove('material_index')


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
    name: str = field(default="coreshell", init=False)
    """ Name of the set """

    def __post_init__(self):
        self.cpp_binding_str: list = [
            'core_diameter',
            'shell_width',
            'core_material_index',
            'shell_material_index',
            'core_index',
            'shell_index',
            'n_medium',
        ]

        self.parameter_dictionnary: dict = dict(
            core_diameter=dict(type=float, value=None),
            shell_width=dict(type=float, value=None),
            core_material=dict(type=complex, value=None),
            core_index=dict(type=complex, value=None),
            shell_material=dict(type=complex, value=None),
            shell_index=dict(type=complex, value=None),
            n_medium=dict(type=float, value=None),
        )

        self.binding_class: type = CppCoreShellSet

        super().__post_init__()

    def asserts_inputs(self) -> None:
        """
        Asserts that core and shell are either defined by index or material

        :returns:   No return
        :rtype:     None
        """
        assert (self.core_material is None) ^ (self.core_index is None), "core_material xor core_index has to be defined"

        assert (self.shell_material is None) ^ (self.shell_index is None), "shell_material xor shell_index has to be defined"

        if self.core_material is not None:
            del self.parameter_dictionnary['core_index']
            self.cpp_binding_str.remove('core_index')
        else:
            del self.parameter_dictionnary['core_material']
            self.cpp_binding_str.remove('core_material_index')

        if self.shell_material is not None:
            del self.parameter_dictionnary['shell_index']
            self.cpp_binding_str.remove('shell_index')
        else:
            del self.parameter_dictionnary['shell_material']
            self.cpp_binding_str.remove('shell_material_index')


@dataclass
class Cylinder(BaseScatterer):
    diameter: Iterable
    """ diameter of the single scatterer in unit of meter. """
    index: Iterable = None
    """ Refractive index of scatterer. """
    material: Iterable = None
    """ Refractive index of scatterer medium. """
    name: str = field(default="cylinder", init=False)
    """ Name of the set """

    def __post_init__(self):
        self.cpp_binding_str: list = [
            'diameter',
            'material_index',
            'index',
            'n_medium',
        ]

        self.parameter_dictionnary: dict = dict(
            diameter=dict(type=float, value=None),
            material=dict(type=complex, value=None),
            index=dict(type=complex, value=None),
            n_medium=dict(type=float, value=None),
        )

        self.binding_class: type = CppCylinderSet

        super().__post_init__()

    def asserts_inputs(self) -> None:
        """
        Asserts that core and shell are either defined by index or material

        :returns:   No return
        :rtype:     None
        """
        assert (self.material is None) ^ (self.index is None), "material xor index has to be defined"

        if self.material is not None:
            del self.parameter_dictionnary['index']
            self.cpp_binding_str.remove('index')
        else:
            del self.parameter_dictionnary['material']
            self.cpp_binding_str.remove('material_index')

# -
