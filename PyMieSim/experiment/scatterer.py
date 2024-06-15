#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from pydantic import ConfigDict
from typing import List, Union, NoReturn, Any

import numpy
from PyMieSim.binary.Sets import CppCoreShellSet, CppCylinderSet, CppSphereSet
from PyMieSim.experiment import measure, parameters
import PyMieSim.experiment.source as source
from PyOptik import Sellmeier, DataMeasurement

config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


@dataclass(config=config_dict)
class BaseScatterer():
    """
    Base class for scatterer objects. This class handles the initialization and setup of
    scatterer parameters for use in PyMieSim simulations.

    """
    mapping = None
    binding_kwargs = None
    binding = None

    def __post_init__(self) -> NoReturn:
        """
        Initializes the scatterer instance by asserting inputs, formatting them, building binding
        arguments, and Units for visualization. This method is automatically called after the
        class has been initialized.
        """
        self.build_binding_kwargs()

    def add_material_index_to_mapping(self, name: str, indexes: numpy.ndarray, materials: numpy.ndarray, data_type: type = object) -> NoReturn:
        """
        Adds material or refractive index details to a mapping dictionary.

        This method is used to create a mapping of material properties to human-readable and accessible formats
        for UI or data outputs. The key in the mapping dictionary is constructed using the provided name.

        Parameters:
            name (str): The base name to use for the keys in the mapping dictionary. This name is used to differentiate between different materials or indices.
            indexes (numpy.ndarray): The array of refractive index values.
            materials (numpy.ndarray): The array of material objects.
            data_type (type): The expected data type of the material or index values. Default is `object`.

        """
        if materials is not None:
            key = f"{name}_material" if name else "material"
            base_values = materials
        else:
            key = f"{name}_index" if name else "index"
            base_values = indexes

        res = getattr(parameters, key)
        res.base_values = res.values = base_values
        self.mapping[key] = res

    def add_material_index_to_binding_kwargs(self, name: str, indexes: numpy.ndarray, materials: numpy.ndarray, data_type: type = object) -> NoReturn:
        """
        Adds either material properties or a refractive index to the binding keyword arguments for the experiment.

        This method validates and processes the material or index information, converting it into a format suitable
        for simulation use, and ensuring that either a material or an index is provided but not both.

        Parameters:
            name (str): The base name for the material or index. This name helps identify the property and is used to handle multiple materials or indices.
            indexes (numpy.ndarray): The array of refractive index values.
            materials (numpy.ndarray): The array of material objects.
            data_type (type): The expected data type of the material or index values. Default is `object`.

        Raises:
            ValueError: If both a material and an index are provided, or if neither is provided.

        Returns:
            NoReturn
        """
        if (materials is not None) == (indexes is not None):
            raise ValueError(f"Either {name} material or {name} index must be provided.")

        if materials:
            key = f"{name}_material" if name else "material"
            wavelength = self.source.wavelength
            materials = numpy.atleast_1d(materials)
            indexes = numpy.asarray([m.get_refractive_index(wavelength) for m in materials])
        else:
            key = f"{name}_index" if name else "index"

        indexes = numpy.atleast_1d(indexes)

        if data_type == float and numpy.any(numpy.iscomplex(indexes)):
            indexes = indexes.real

        self.binding_kwargs[key] = indexes.astype(data_type)


@dataclass(config=config_dict)
class Sphere(BaseScatterer):
    """
    A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

    This class provides specific implementations for setting up and binding spherical scatterers
    with their properties to a simulation environment. It extends the `BaseScatterer` class by
    adding spherical-specific attributes and methods for handling simulation setups.

    Attributes:
        source (Union[experiment.source.Gaussian, experiment.source.PlaneWave]): Light source configuration for the simulation.
        diameter (List): Diameter(s) of the spherical scatterers in meters.
        medium_index (List, optional): Refractive index or indices of the medium surrounding the scatterers.
        medium_material (List, optional): Material(s) defining the medium, used if `medium_index` is not provided.
        index (List, optional): Refractive index or indices of the spherical scatterers themselves.
        material (List, optional): Material(s) of the scatterers, used if `index` is not provided.
        name (str): Name identifier for the scatterer type, defaulted to 'sphere' and not intended for initialization.
    """
    source: Union[source.Gaussian, source.PlaneWave]
    diameter: Union[numpy.ndarray, List[float], float]
    medium_index: Union[numpy.ndarray, List[float], float, None] = None
    medium_material: Union[List[Sellmeier | DataMeasurement], Sellmeier | DataMeasurement, None] = None
    index: Union[numpy.ndarray, List[Any], Any] = None
    material: Union[List[Sellmeier | DataMeasurement], Sellmeier | DataMeasurement, None] = None

    available_measure_list = measure.__sphere__

    def __post_init__(self):
        """
        Extends the base class post-initialization process by setting up additional properties specific to spherical scatterers.
        Initializes a mapping dictionary to support visualization and other operations.
        """
        self.mapping = dict.fromkeys(['diameter'])

        self.diameter = numpy.atleast_1d(self.diameter).astype(float)

        super(Sphere, self).__post_init__()

    def build_binding_kwargs(self) -> NoReturn:
        """
        Constructs the keyword arguments necessary for the C++ binding interface, specifically tailored for spherical scatterers.
        This includes processing material indices and organizing them into a structured dictionary for simulation interaction.

        This method automatically configures the `binding_kwargs` attribute with appropriately formatted values.
        """
        self.binding_kwargs = dict(diameter=self.diameter)

        self.add_material_index_to_binding_kwargs(
            name=None,
            indexes=self.index,
            materials=self.material,
            data_type=complex
        )

        self.add_material_index_to_binding_kwargs(
            name='medium',
            indexes=self.medium_index,
            materials=self.medium_material,
            data_type=float
        )

        self.binding = CppSphereSet(**self.binding_kwargs)

    def get_datavisual_table(self) -> NoReturn:
        """
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns:
            list: A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """
        self.mapping['diameter'] = parameters.diameter
        self.mapping['diameter'].set_base_values(self.diameter)

        self.add_material_index_to_mapping(
            name=None,
            indexes=self.index,
            materials=self.material,
        )
        self.add_material_index_to_mapping(
            name='medium',
            indexes=self.medium_index,
            materials=self.medium_material,
        )

        return [v for k, v in self.mapping.items() if v is not None]


@dataclass(config=config_dict)
class CoreShell(BaseScatterer):
    """
    A data class representing a core-shell scatterer configuration used in PyMieSim simulations.

    This class facilitates the setup and manipulation of core-shell scatterers by providing structured
    attributes and methods that ensure the scatterers are configured correctly for simulations.
    It extends the BaseScatterer class, adding specific attributes and methods relevant to core-shell geometries.

    Attributes:
        source (Union[experiment.source.Gaussian, experiment.source.PlaneWave]): Light source configuration for the simulation.
        core_diameter (Union[List[float], float]): Diameters of the core components in meters.
        shell_width (Union[List[float], float]): Thicknesses of the shell components in meters.
        medium_index (List, optional): Refractive index or indices of the medium where the scatterers are placed.
        medium_material (List, optional): Material(s) defining the medium, used if `medium_index` is not provided.
        core_index (List, optional): Refractive index or indices of the core.
        shell_index (List, optional): Refractive index or indices of the shell.
        core_material (List, optional): Material(s) of the core, used if `core_index` is not provided.
        shell_material (List, optional): Material(s) of the shell, used if `shell_index` is not provided.
        name (str): An identifier for the scatterer type, defaulted to 'coreshell' and not intended for initialization.
    """
    source: Union[source.Gaussian, source.PlaneWave]
    core_diameter: Union[numpy.ndarray, List[float], float]
    shell_width: Union[numpy.ndarray, List[float], float]
    medium_index: Union[numpy.ndarray, List[float], float, None] = None
    medium_material: Union[List[Sellmeier | DataMeasurement], Sellmeier | DataMeasurement, None] = None
    shell_index: Union[numpy.ndarray, List[Any], Any, None] = None
    core_material: Union[List[Sellmeier | DataMeasurement], Sellmeier | DataMeasurement, None] = None
    core_index: Union[numpy.ndarray, List[Any], Any, None] = None
    shell_material: Union[List[Sellmeier | DataMeasurement], Sellmeier | DataMeasurement, None] = None

    available_measure_list = measure.__coreshell__

    def __post_init__(self):
        """
        Extends the BaseScatterer post-initialization by setting up additional mappings specific to core-shell scatterers.
        Initializes mappings for visualizing and interacting with scatterer properties.
        """
        self.mapping = dict.fromkeys([
            'core_diameter',
            'shell_width',
        ])

        self.core_diameter = numpy.atleast_1d(self.core_diameter).astype(float)
        self.shell_width = numpy.atleast_1d(self.shell_width).astype(float)

        super(CoreShell, self).__post_init__()

    def build_binding_kwargs(self) -> NoReturn:
        """
        Assembles the keyword arguments necessary for C++ binding, tailored for core-shell scatterers.
        Prepares structured data from scatterer properties for efficient simulation integration.

        This function populates `binding_kwargs` with values formatted appropriately for the C++ backend used in simulations.
        """
        self.binding_kwargs = dict(core_diameter=self.core_diameter, shell_width=self.shell_width)

        self.add_material_index_to_binding_kwargs(
            name='core',
            indexes=self.core_index,
            materials=self.core_material,
            data_type=complex
        )

        self.add_material_index_to_binding_kwargs(
            name='shell',
            indexes=self.shell_index,
            materials=self.shell_material,
            data_type=complex
        )

        self.add_material_index_to_binding_kwargs(
            name='medium',
            indexes=self.medium_index,
            materials=self.medium_material,
            data_type=float
        )

        self.binding = CppCoreShellSet(**self.binding_kwargs)

    def get_datavisual_table(self) -> NoReturn:
        """
        Generates a list of data visualizations for the scatterer's properties, which can be used in user interfaces or reports.
        Each property is formatted into a user-friendly structure, making it easier to visualize and understand.

        Returns:
            list: A collection of formatted representations for the scatterer properties.
        """
        self.mapping['core_diameter'] = parameters.core_diameter
        self.mapping['core_diameter'].set_base_values(self.binding_kwargs.get('core_diameter'))

        self.mapping['shell_width'] = parameters.shell_width
        self.mapping['shell_width'].set_base_values(self.binding_kwargs.get('shell_width'))

        self.add_material_index_to_mapping(
            name='core',
            indexes=self.core_index,
            materials=self.core_material,
        )

        self.add_material_index_to_mapping(
            name='shell',
            indexes=self.shell_index,
            materials=self.shell_material,
        )

        self.add_material_index_to_mapping(
            name='medium',
            indexes=self.medium_index,
            materials=self.medium_material,
        )

        return [v for k, v in self.mapping.items() if v is not None]


@dataclass(config=config_dict)
class Cylinder(BaseScatterer):
    """
    Represents a cylindrical scatterer configuration for PyMieSim simulations.

    Attributes:
        source (Union[experiment.source.Gaussian, experiment.source.PlaneWave]): Light source configuration for the simulation.
        diameter (List): Diameter(s) of the cylinder in meters.
        height (List): Height(s) of the cylinder in meters.
        index (List, optional): Refractive index of the cylinder.
        material (List, optional): Material(s) of the cylinder, used if `index` is not provided.
    """
    source: Union[source.Gaussian, source.PlaneWave]
    diameter: Union[numpy.ndarray, List[float], float]
    medium_index: Union[numpy.ndarray, List[float], float, None] = None
    medium_material: Union[List[Sellmeier | DataMeasurement], Sellmeier | DataMeasurement, None] = None
    index: Union[numpy.ndarray, List[Any], Any, None] = None
    material: Union[List[Sellmeier | DataMeasurement], Sellmeier | DataMeasurement, None] = None

    available_measure_list = measure.__cylinder__

    def __post_init__(self):
        self.mapping = dict.fromkeys(['diameter'])
        self.diameter = numpy.atleast_1d(self.diameter).astype(float)

        super(Cylinder, self).__post_init__()

    def build_binding_kwargs(self) -> NoReturn:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns:
            None
        """
        self.binding_kwargs = dict(diameter=self.diameter)

        self.add_material_index_to_binding_kwargs(
            name=None,
            indexes=self.index,
            materials=self.material,
            data_type=complex
        )

        self.add_material_index_to_binding_kwargs(
            name='medium',
            indexes=self.medium_index,
            materials=self.medium_material,
            data_type=float
        )

        self.binding = CppCylinderSet(**self.binding_kwargs)

    def get_datavisual_table(self) -> List:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['diameter'] = parameters.diameter
        self.mapping['diameter'].set_base_values(self.binding_kwargs.get('diameter'))

        self.add_material_index_to_mapping(
            name=None,
            indexes=self.index,
            materials=self.material,
        )
        self.add_material_index_to_mapping(
            name='medium',
            indexes=self.medium_index,
            materials=self.medium_material,
        )

        return [v for k, v in self.mapping.items() if v is not None]

# -
