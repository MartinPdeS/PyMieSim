#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from dataclasses import dataclass

from DataVisual import Xparameter

import PyMieSim.datavisual_x_parameters as Kwargs

from PyMieSim.binary.Sets import CppCoreShellSet, CppCylinderSet, CppSphereSet


@dataclass
class CoreShell():
    core_diameter: list
    """ diameter of the core of the single scatterer [m]. """
    shell_width: list
    """ Width of the shell of the single scatterer [m]. """
    core_index: tuple = tuple()
    """ Refractive index of the core of the scatterer. """
    shell_index: tuple = tuple()
    """ Refractive index of the shell of the scatterer. """
    core_material: tuple = tuple()
    """ Core material of which the scatterer is made of. Only if core_index is not specified.  """
    shell_material: tuple = tuple()
    """ Shell material of which the scatterer is made of. Only if shell_index is not specified.  """
    n_medium: list = 1.0
    """ Refractive index of scatterer medium. """
    name: str = 'coreshell'
    """name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        super().__init__()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.n_medium = Xparameter(values=self.n_medium, **Kwargs.n_medium)

        self.bounded_core = True if len(self.core_material) != 0 else False

        self.bounded_shell = True if len(self.shell_material) != 0 else False

        self.core_diameter = Xparameter(values=self.core_diameter, **Kwargs.core_diameter)

        self.shell_width = Xparameter(values=self.shell_width, **Kwargs.shell_width)

        self.core_material = Xparameter(values=self.core_material, **Kwargs.core_material)

        self.shell_material = Xparameter(values=self.shell_material, **Kwargs.shell_material)

        self.core_index = Xparameter(values=self.core_index, **Kwargs.core_index)

        self.shell_index = Xparameter(values=self.shell_index, **Kwargs.shell_index)

    def bind_to_experiment(self, experiment):
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_coreshell(self.binding)

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

    def bind_core_material_shell_material(self, source) -> None:
        core_material = [
            material.GetRI(source.values) for material in self.core_material
        ]

        shell_material = [
            material.GetRI(source.values) for material in self.shell_material
        ]

        core_material = numpy.asarray(core_material).astype(complex)
        shell_material = numpy.asarray(shell_material).astype(complex)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_material=core_material,
            shell_material=shell_material,
            n_medium=self.n_medium.values,
        )

    def bind_core_material_shell_index(self, source) -> None:
        core_material = [
            material.GetRI(source.values) for material in self.core_material
        ]

        core_material = numpy.asarray(core_material).astype(complex)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_material=core_material,
            shell_index=self.shell_index.values,
            n_medium=self.n_medium.values,
        )

    def bind_index_material_shell_material(self, source) -> None:
        shell_material = [
            material.GetRI(source.values) for material in self.shell_material
        ]

        shell_material = numpy.asarray(shell_material).astype(complex)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_index=self.core_index.values,
            shell_material=shell_material,
            n_medium=self.n_medium.values,
        )

    def bind_index_material_shell_index(self, source) -> None:
        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_index=self.core_index.values,
            shell_index=self.shell_index.values,
            n_medium=self.n_medium.values,
        )

    def evaluate_index_material(self, source) -> None:
        if self.bounded_core and self.bounded_shell:
            self.bind_core_material_shell_material(source)

        if self.bounded_core and not self.bounded_shell:
            self.bind_core_material_shell_index(source)

        if not self.bounded_core and self.bounded_shell:
            self.bind_index_material_shell_material(source)

        if not self.bounded_core and not self.bounded_shell:
            self.bind_index_material_shell_index(source)

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.core_diameter = numpy.atleast_1d(self.core_diameter).astype(float)
        self.shell_width = numpy.atleast_1d(self.shell_width).astype(float)

        self.core_index = numpy.atleast_1d(self.core_index).astype(numpy.complex128)
        self.shell_index = numpy.atleast_1d(self.shell_index).astype(numpy.complex128)

        self.n_medium = numpy.atleast_1d(self.n_medium).astype(float)

        self.core_material = numpy.atleast_1d(self.core_material)
        self.shell_material = numpy.atleast_1d(self.shell_material)

    def append_to_table(self, table) -> list:
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
class Sphere():
    diameter: list
    """ diameter of the single scatterer in unit of meter. """
    index: list = tuple()
    """ Refractive index of scatterer. """
    material: list = tuple()
    """ material of which the scatterer is made of. Only if index is not specified. """
    n_medium: list = 1.0
    """ Refractive index of scatterer medium. """
    name: str = 'sphere'
    """name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        super().__init__()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.n_medium = Xparameter(values=self.n_medium, **Kwargs.n_medium)

        self.bounded_index = True if len(self.material) != 0 else False

        self.diameter = Xparameter(values=self.diameter, **Kwargs.diameter)

        self.material = Xparameter(values=self.material, **Kwargs.material)

        self.index = Xparameter(values=self.index, **Kwargs.index)

    def bind_to_experiment(self, experiment) -> None:
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_sphere(self.binding)

    def evaluate_index_material(self, Source) -> None:
        if self.bounded_index:
            material = [
                material.GetRI(Source.values) for material in self.material
            ]

            material_index = numpy.asarray(material).astype(complex)

            self.binding = CppSphereSet(
                diameter=self.diameter.values.astype(float),
                material=material_index,
                n_medium=self.n_medium.values.astype(float)
            )

        else:
            self.binding = CppSphereSet(
                diameter=self.diameter.values,
                index=self.index.values,
                n_medium=self.n_medium.values
            )

    def append_to_table(self, table) -> list:
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

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.diameter = numpy.atleast_1d(self.diameter).astype(float)

        self.index = numpy.atleast_1d(self.index).astype(numpy.complex128)

        self.n_medium = numpy.atleast_1d(self.n_medium).astype(float)

        self.material = numpy.atleast_1d(self.material)


@dataclass
class Cylinder():
    diameter: list
    """ diameter of the single scatterer in unit of meter. """
    index: tuple = tuple()
    """ Refractive index of scatterer. """
    material: tuple = tuple()
    """ Refractive index of scatterer medium. """
    n_medium: list = 1.0
    """ material of which the scatterer is made of. Only if index is not specified. """
    name: str = 'cylinder'
    """name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        super().__init__()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.n_medium = Xparameter(values=self.n_medium, **Kwargs.n_medium)

        self.bounded_index = True if len(self.material) != 0 else False

        self.diameter = Xparameter(values=self.diameter, **Kwargs.diameter)

        self.material = Xparameter(values=self.material, **Kwargs.material)

        self.index = Xparameter(values=self.index, **Kwargs.index)

    def bind_to_experiment(self, experiment):
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_cylinder(self.binding)

    def evaluate_index_material(self, source):
        if self.bounded_index:
            material = numpy.asarray([material.GetRI(source.values) for material in self.material])

            self.binding = CppCylinderSet(
                diameter=self.diameter.values,
                material=material,
                n_medium=self.n_medium.values
            )

        else:
            self.binding = CppCylinderSet(
                diameter=self.diameter.values,
                index=self.index.values,
                n_medium=self.n_medium.values
            )

    def append_to_table(self, table):
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

    def format_inputs(self):
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.diameter = numpy.atleast_1d(self.diameter).astype(float)

        self.index = numpy.atleast_1d(self.index).astype(numpy.complex128)

        self.n_medium = numpy.atleast_1d(self.n_medium).astype(float)

        self.material = numpy.atleast_1d(self.material)

# -
