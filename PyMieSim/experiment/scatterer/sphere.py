#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from pydantic import field_validator
from typing import List, Union, Any, Optional

import numpy
from PyMieSim.binary.SetsInterface import CppSphereSet
from PyMieSim.experiment import measure
import PyMieSim.experiment.source as source
from PyOptik.base_class import BaseMaterial
from PyMieSim.units import Quantity, meter
from PyMieSim.experiment.scatterer.base_class import BaseScatterer, config_dict

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
        medium_material (List, optional): BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
        index (List, optional): Refractive index or indices of the spherical scatterers themselves.
        material (List, optional): BaseMaterial(s) of the scatterers, used if `index` is not provided.
        name (str): Name identifier for the scatterer type, defaulted to 'sphere' and not intended for initialization.
    """
    source: Union[source.Gaussian, source.PlaneWave]
    diameter: Quantity
    medium_index: Optional[Quantity] = None
    medium_material: Optional[List[BaseMaterial] | BaseMaterial] = None
    index: Optional[Quantity] = None
    material: Optional[List[BaseMaterial] | BaseMaterial] = None

    available_measure_list = measure.__sphere__

    @field_validator('diameter', 'medium_index', 'medium_material', 'index', 'material', mode='before')
    def validate_array(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    @field_validator('diameter', mode='before')
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
        Constructs the keyword arguments necessary for the C++ binding interface, specifically tailored for spherical scatterers.
        This includes processing material indices and organizing them into a structured dictionary for simulation interaction.

        This method automatically configures the `binding_kwargs` attribute with appropriately formatted values.
        """
        self.mapping = {}

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

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppSphereSet(**binding_kwargs)