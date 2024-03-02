#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from typing import NoReturn
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
    """
    Base class for scatterer objects. This class handles the initialization and setup of
    scatterer parameters for use in PyMieSim simulations.

    Attributes:
        n_medium (Iterable): Refractive index of the medium in which the scatterers are placed.
        source_set (Union[Gaussian, PlaneWave]): Light source configuration for the simulation.
    """
    n_medium: Iterable
    source_set: Gaussian | PlaneWave

    def __post_init__(self) -> NoReturn:
        """
        Initializes the scatterer instance by asserting inputs, formatting them, building binding
        arguments, and Xparameters for visualization. This method is automatically called after the
        class has been initialized.

        Returns:
            NoReturn
        """
        self.validate_material_or_index()

        self.format_inputs()

        self.build_binding_kwargs()

        self.build_x_parameters()

    def _validate_and_cleanup(self, part: str = '') -> NoReturn:
        """
        Validates and cleans up parameters for either the core or the shell.

        Args:
            part (str): Specifies the part to validate, either 'core' or 'shell' or ''.
        """
        material = getattr(self, f"{part}material")
        index = getattr(self, f"{part}index")

        if material is not None and index is not None:
            raise ValueError(f"Either {part} material or index must be provided, not both.")
        if material is None and index is None:
            raise ValueError(f"One of {part} material or index must be provided.")

        # Cleanup parameters
        if material is not None:
            del self.parameter_dictionnary[f'{part}index']
            self.cpp_binding_str.remove(f'{part}index')
        else:
            del self.parameter_dictionnary[f'{part}material']
            self.cpp_binding_str.remove(f'{part}material_index')

    def format_inputs(self) -> NoReturn:
        """
        Formats the input attributes into numpy arrays for further processing and ensures compatibility
        with C++ binding. This method standardizes the input parameters for simulation.

        Returns:
            None
        """
        for parameter_str in self.parameter_dictionnary.keys():
            parameter = getattr(self, parameter_str)

            parameter = numpy.atleast_1d(parameter)

            setattr(self, parameter_str, parameter)

    def bind_to_experiment(self, experiment: Setup) -> NoReturn:
        """
        Binds the scatterer to a specific experiment setup, enabling its properties to be evaluated within
        the given experimental context.

        Parameters:
            experiment (Setup): The experiment setup to which the scatterer will be bound.

        Returns:
            None
        """
        method_str = 'set_' + self.name

        getattr(experiment.binding, method_str)(self.binding)

    def build_binding_kwargs(self) -> NoReturn:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns:
            None
        """
        self.evaluate_index_material()

        self.binding_kwargs = dict()

        for parameter_str in self.cpp_binding_str:
            values = getattr(self, parameter_str)

            self.binding_kwargs[parameter_str] = values

    def build_x_parameters(self) -> NoReturn:
        """
        Constructs Xparameters for inclusion in the XTable for DataVisual, facilitating the visualization
        of the scatterer's properties.

        Returns:
            None
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
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        return [*table, *self.x_table]

    def initialize_binding(self) -> NoReturn:
        """
        Initializes the C++ binding for the scatterer, setting up the interface for simulation within the
        C++ computational framework.

        Returns:
            None
        """
        self.binding = self.binding_class(**self.binding_kwargs)

    def evaluate_index_material(self) -> NoReturn:
        """
        Evaluates the indices of materials for all scatterer sets that are specified by their material
        properties. This is crucial for accurately simulating the scatterer's interaction with light.

        Returns:
            None
        """
        for parameter_str, value in self.parameter_dictionnary.items():
            if parameter_str.endswith('material'):
                materials = getattr(self, parameter_str)

                material_index = [
                    material.get_refractive_index(self.source_set.wavelength.values) for material in materials
                ]

                material_index = numpy.asarray(material_index).astype(complex)

                setattr(self, parameter_str + '_index', material_index)


@dataclass
class Sphere(BaseScatterer):
    """
    Represents a spherical scatterer configuration for PyMieSim simulations.

    Attributes:
        diameter (Iterable): Diameter(s) of the scatterers in meters.
        index (Iterable, optional): Refractive index of the scatterers.
        material (Iterable, optional): Material(s) of the scatterers, used if `index` is not provided.
    """
    diameter: Iterable
    index: Iterable | None = None
    material: Iterable | None = None

    name: str = field(default="sphere", init=False)

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

    def validate_material_or_index(self) -> NoReturn:
        """
        Validates the inputs for the CoreShell scatterer, ensuring that both core and shell are defined
        either by their material or by their refractive index, but not both simultaneously. This ensures
        accurate modeling of core-shell scatterers.

        Returns:
            NoReturn
        """
        self._validate_and_cleanup()


@dataclass
class CoreShell(BaseScatterer):
    """
    Represents a core-shell scatterer configuration for PyMieSim simulations.

    Attributes:
        core_diameter (Iterable): Diameter(s) of the core in meters.
        shell_thickness (Iterable): Thickness(es) of the shell in meters.
        core_index (Iterable, optional): Refractive index of the core.
        shell_index (Iterable, optional): Refractive index of the shell.
        core_material (Iterable, optional): Material(s) of the core, used if `core_index` is not provided.
        shell_material (Iterable, optional): Material(s) of the shell, used if `shell_index` is not provided.
    """
    core_diameter: Iterable
    shell_width: Iterable
    core_index: Iterable | None = None
    shell_index: Iterable | None = None
    core_material: Iterable | None = None
    shell_material: Iterable | None = None

    name: str = field(default="coreshell", init=False)

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

    def validate_material_or_index(self) -> NoReturn:
        """
        Validates the inputs for the CoreShell scatterer, ensuring that both core and shell are defined
        either by their material or by their refractive index, but not both simultaneously. This ensures
        accurate modeling of core-shell scatterers.

        Returns:
            NoReturn
        """
        self._validate_and_cleanup('core_')
        self._validate_and_cleanup('shell_')


@dataclass
class Cylinder(BaseScatterer):
    """
    Represents a cylindrical scatterer configuration for PyMieSim simulations.

    Attributes:
        diameter (Iterable): Diameter(s) of the cylinder in meters.
        height (Iterable): Height(s) of the cylinder in meters.
        index (Iterable, optional): Refractive index of the cylinder.
        material (Iterable, optional): Material(s) of the cylinder, used if `index` is not provided.
    """
    diameter: Iterable
    index: Iterable | None = None
    material: Iterable | None = None

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

    def validate_material_or_index(self) -> NoReturn:
        """
        Validates the inputs for the CoreShell scatterer, ensuring that both core and shell are defined
        either by their material or by their refractive index, but not both simultaneously. This ensures
        accurate modeling of core-shell scatterers.

        Returns:
            NoReturn
        """
        self._validate_and_cleanup()

# -
