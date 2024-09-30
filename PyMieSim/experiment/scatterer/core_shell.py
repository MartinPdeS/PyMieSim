#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from pydantic import field_validator
from typing import List, Union, Any, Optional

import numpy
from PyMieSim.binary.SetsInterface import CppCoreShellSet
from PyMieSim.experiment import measure
import PyMieSim.experiment.source as source
from PyOptik.base_class import BaseMaterial
from PyMieSim.units import Quantity, meter
from PyMieSim.experiment.scatterer.base_class import BaseScatterer, config_dict


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
        medium_material (List, optional): BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
        core_index (List, optional): Refractive index or indices of the core.
        shell_index (List, optional): Refractive index or indices of the shell.
        core_material (List, optional): BaseMaterial(s) of the core, used if `core_index` is not provided.
        shell_material (List, optional): BaseMaterial(s) of the shell, used if `shell_index` is not provided.
        name (str): An identifier for the scatterer type, defaulted to 'coreshell' and not intended for initialization.
    """
    source: Union[source.Gaussian, source.PlaneWave]
    core_diameter: Quantity
    shell_width: Quantity
    medium_index: Optional[Quantity] = None
    medium_material: Optional[List[BaseMaterial] | BaseMaterial] = None
    shell_index: Optional[Quantity] = None
    core_material: Optional[List[BaseMaterial] | BaseMaterial] = None
    core_index: Optional[Quantity] = None
    shell_material: Optional[List[BaseMaterial] | BaseMaterial] = None

    available_measure_list = measure.__coreshell__

    @field_validator('core_diameter', 'shell_width', 'core_index', 'shell_index', 'core_material', 'shell_material', 'medium_index', 'medium_material', mode='before')
    def validate_array(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    @field_validator('core_diameter', 'shell_wdith', mode='before')
    def validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return numpy.atleast_1d(value)

    def __post_init__(self) -> None:
        """
        Assembles the keyword arguments necessary for C++ binding, tailored for core-shell scatterers.
        Prepares structured data from scatterer properties for efficient simulation integration.

        This function populates `binding_kwargs` with values formatted appropriately for the C++ backend used in simulations.
        """
        self.mapping = {}

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

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppCoreShellSet(**binding_kwargs)
