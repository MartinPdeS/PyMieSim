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

from DataVisual import units
from PyMieSim.binary.Sets import CppCoreShellSet, CppCylinderSet, CppSphereSet
from PyMieSim import measure


@dataclass
class BaseScatterer():
    """
    Base class for scatterer objects. This class handles the initialization and setup of
    scatterer parameters for use in PyMieSim simulations.

    Attributes:
        medium_index (Iterable): Refractive index of the medium in which the scatterers are placed.
        source_set (Union[Gaussian, PlaneWave]): Light source configuration for the simulation.
    """
    source_set: Gaussian | PlaneWave

    def __post_init__(self) -> NoReturn:
        """
        Initializes the scatterer instance by asserting inputs, formatting them, building binding
        arguments, and Units for visualization. This method is automatically called after the
        class has been initialized.

        Returns:
            NoReturn
        """
        self.validate_material_or_index()

        self.build_binding_kwargs()

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

    def add_material_index_to_mapping(self, name: str = None) -> NoReturn:
        detached_material_name = f"{name} material" if name else "material"
        attached_material_name = detached_material_name.replace(' ', '_').lower()

        detached_index_name = f"{name} index" if name else "index"
        attached_index_name = detached_index_name.replace(' ', '_').lower()

        if getattr(self, attached_material_name):
            self.mapping[attached_material_name] = units.Custom(
                long_label=detached_material_name,
                short_label=attached_material_name,
                value_representation=getattr(self, attached_material_name),
                base_values=getattr(self, attached_material_name),
                use_prefix=False,
            )

        else:
            self.mapping[attached_index_name] = units.Index(
                long_label=detached_index_name,
                short_label=attached_index_name,
                base_values=self.binding_kwargs.get(attached_index_name),
                string_format='.2f'
            )


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
    medium_index: Iterable | None = None
    medium_material: Iterable | None = None
    index: Iterable | None = None
    material: Iterable | None = None
    name: str = field(default="sphere", init=False)
    available_measure_list = measure.__sphere__

    def __post_init__(self):
        self.mapping: dict = {
            'diameter': None,
            'index': None,
            'material': None,
            'medium_index': None,
            'medium_material': None
        }

        super().__post_init__()

    def build_binding_kwargs(self) -> NoReturn:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns:
            None
        """
        self.binding_kwargs = dict(
            diameter=numpy.atleast_1d(self.diameter).astype(float),
        )

        if self.material:
            self.binding_kwargs['material'] = numpy.asarray([
                mat.get_refractive_index(self.source_set.wavelength) for mat in numpy.atleast_1d(self.material)
            ]).astype(complex)
        else:
            self.binding_kwargs['index'] = numpy.atleast_1d(self.index).astype(complex)

        if self.medium_material:
            self.binding_kwargs['medium_material'] = numpy.asarray([
                mat.get_refractive_index(self.source_set.wavelength) for mat in numpy.atleast_1d(self.medium_material)
            ]).real.astype(float)
        else:
            self.binding_kwargs['medium_index'] = numpy.atleast_1d(self.medium_index).astype(float)

        self.binding = CppSphereSet(**self.binding_kwargs)

    def get_datavisual_table(self) -> NoReturn:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['diameter'] = units.Length(
            long_label='Scatterer diameter',
            short_label='diameter',
            base_values=self.binding_kwargs.get('diameter'),
            string_format='.1f'
        )

        self.add_material_index_to_mapping(name=None)
        self.add_material_index_to_mapping(name='Medium')

        return [v for k, v in self.mapping.items() if v is not None]

    def validate_material_or_index(self) -> NoReturn:
        """
        Validates the inputs for the CoreShell scatterer, ensuring that both core and shell are defined
        either by their material or by their refractive index, but not both simultaneously. This ensures
        accurate modeling of core-shell scatterers.

        Returns:
            NoReturn
        """
        if self.material is not None and self.index is not None:
            raise ValueError(f"Either material or index must be provided, not both.")
        if self.material is None and self.index is None:
            raise ValueError(f"One of material or index must be provided.")

        if self.medium_material is not None and self.medium_index is not None:
            raise ValueError(f"Either material or index must be provided, not both.")
        if self.medium_material is None and self.medium_index is None:
            raise ValueError(f"One of material or index must be provided.")


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
    medium_index: Iterable | None = None
    medium_material: Iterable | None = None
    core_index: Iterable | None = None
    shell_index: Iterable | None = None
    core_material: Iterable | None = None
    shell_material: Iterable | None = None
    name: str = field(default="coreshell", init=False)

    available_measure_list = measure.__coreshell__

    def __post_init__(self):
        self.mapping = {
            'core_diameter': None,
            'shell_diameter': None,
            'core_index': None,
            'core_material': None,
            'shell_index': None,
            'shell_material': None,
            'medium_index': None,
            'medium_material': None
        }

        self.binding_class: type = CppCoreShellSet

        super().__post_init__()

    def build_binding_kwargs(self) -> NoReturn:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns:
            None
        """
        self.binding_kwargs = dict(
            core_diameter=numpy.atleast_1d(self.core_diameter).astype(float),
            shell_width=numpy.atleast_1d(self.shell_width).astype(float),
        )

        if self.core_material:
            self.binding_kwargs['core_material'] = numpy.asarray([
                mat.get_refractive_index(self.source_set.wavelength) for mat in numpy.atleast_1d(self.core_material)
            ]).astype(complex)

        else:
            self.binding_kwargs['core_index'] = numpy.atleast_1d(self.core_index).astype(complex)

        if self.shell_material:
            self.binding_kwargs['shell_material'] = numpy.asarray([
                mat.get_refractive_index(self.source_set.wavelength) for mat in numpy.atleast_1d(self.shell_material)
            ]).astype(complex)
        else:
            self.binding_kwargs['shell_index'] = numpy.atleast_1d(self.shell_index).astype(complex)

        if self.medium_material:
            self.binding_kwargs['medium_material'] = numpy.asarray([
                mat.get_refractive_index(self.source_set.wavelength) for mat in numpy.atleast_1d(self.medium_material)
            ]).real.astype(float)
        else:
            self.binding_kwargs['medium_index'] = numpy.atleast_1d(self.medium_index).astype(float)

        self.binding = CppCoreShellSet(**self.binding_kwargs)

    def get_datavisual_table(self) -> NoReturn:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['core_diameter'] = units.Length(
            long_label='Core diameter',
            short_label='core_diameter',
            base_values=self.binding_kwargs.get('core_diameter'),
            string_format='.2f'
        )

        self.mapping['shell_width'] = units.Length(
            long_label='Shell width',
            short_label='shell_width',
            base_values=self.binding_kwargs.get('shell_width'),
            string_format='.2f'
        )

        self.add_material_index_to_mapping(name='Core')
        self.add_material_index_to_mapping(name='Shell')
        self.add_material_index_to_mapping(name='Medium')

        return [v for k, v in self.mapping.items() if v is not None]

    def validate_material_or_index(self) -> NoReturn:
        """
        Validates the inputs for the CoreShell scatterer, ensuring that both core and shell are defined
        either by their material or by their refractive index, but not both simultaneously. This ensures
        accurate modeling of core-shell scatterers.

        Returns:
            NoReturn
        """
        if self.core_material is not None and self.core_index is not None:
            raise ValueError(f"Either core material or index must be provided, not both.")
        if self.core_material is None and self.core_index is None:
            raise ValueError(f"One of core material or index must be provided.")

        if self.shell_material is not None and self.shell_index is not None:
            raise ValueError(f"Either shell material or index must be provided, not both.")
        if self.shell_material is None and self.shell_index is None:
            raise ValueError(f"One of shell material or index must be provided.")

        if self.medium_material is not None and self.medium_index is not None:
            raise ValueError(f"Either material or index must be provided, not both.")
        if self.medium_material is None and self.medium_index is None:
            raise ValueError(f"One of material or index must be provided.")


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
    medium_index: Iterable | None = None
    medium_material: Iterable | None = None
    index: Iterable | None = None
    material: Iterable | None = None
    name: str = field(default="cylinder", init=False)

    available_measure_list = measure.__cylinder__

    def __post_init__(self):
        self.mapping = {
            'diameter': None,
            'index': None,
            'material': None,
            'medium_index': None,
            'medium_material': None
        }

        super().__post_init__()

    def build_binding_kwargs(self) -> NoReturn:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns:
            None
        """
        self.binding_kwargs = dict(
            diameter=numpy.atleast_1d(self.diameter).astype(float),
        )

        if self.material:
            self.binding_kwargs['material'] = numpy.asarray([
                mat.get_refractive_index(self.source_set.wavelength) for mat in numpy.atleast_1d(self.material)
            ]).astype(complex)
        else:
            self.binding_kwargs['index'] = numpy.atleast_1d(self.index).astype(complex)

        if self.medium_material:
            self.binding_kwargs['medium_material'] = numpy.asarray([
                mat.get_refractive_index(self.source_set.wavelength) for mat in numpy.atleast_1d(self.medium_material)
            ]).real.astype(float)
        else:
            self.binding_kwargs['medium_index'] = numpy.atleast_1d(self.medium_index).astype(float)

        self.binding = CppCylinderSet(**self.binding_kwargs)

    def get_datavisual_table(self) -> NoReturn:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['diameter'] = units.Length(
            long_label='Scatterer diameter',
            short_label='diameter',
            base_values=self.binding_kwargs.get('diameter'),
            string_format='.1f'
        )

        self.add_material_index_to_mapping(name=None)
        self.add_material_index_to_mapping(name='Medium')

        return [v for k, v in self.mapping.items() if v is not None]

    def validate_material_or_index(self) -> NoReturn:
        """
        Validates the inputs for the CoreShell scatterer, ensuring that both core and shell are defined
        either by their material or by their refractive index, but not both simultaneously. This ensures
        accurate modeling of core-shell scatterers.

        Returns:
            NoReturn
        """
        if self.material is not None and self.index is not None:
            raise ValueError(f"Either material or index must be provided, not both.")
        if self.material is None and self.index is None:
            raise ValueError(f"One of material or index must be provided.")

        if self.medium_material is not None and self.medium_index is not None:
            raise ValueError(f"Either material or index must be provided, not both.")
        if self.medium_material is None and self.medium_index is None:
            raise ValueError(f"One of material or index must be provided.")
# -
