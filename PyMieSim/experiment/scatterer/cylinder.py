#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from pydantic import validator
from typing import List, Union, Any, Optional
import pint_pandas

import numpy
from PyMieSim.binary.SetsInterface import CppCylinderSet
from PyMieSim.experiment import measure, parameters
import PyMieSim.experiment.source as source
from PyOptik.base_class import BaseMaterial
from PyMieSim.units import Quantity, meter
from PyMieSim.experiment.scatterer.base_class import BaseScatterer, config_dict


@dataclass(config=config_dict)
class Cylinder(BaseScatterer):
    """
    Represents a cylindrical scatterer configuration for PyMieSim simulations.

    Attributes:
        source (Union[experiment.source.Gaussian, experiment.source.PlaneWave]): Light source configuration for the simulation.
        diameter (List): Diameter(s) of the cylinder in meters.
        height (List): Height(s) of the cylinder in meters.
        index (List, optional): Refractive index of the cylinder.
        material (List, optional): BaseMaterial(s) of the cylinder, used if `index` is not provided.
    """
    source: Union[source.Gaussian, source.PlaneWave]
    diameter: Quantity
    medium_index: Optional[Quantity] = None
    medium_material: Optional[List[BaseMaterial] | BaseMaterial] = None
    index: Optional[Quantity] = None
    material: Optional[List[BaseMaterial] | BaseMaterial] = None

    available_measure_list = measure.__cylinder__

    @validator('diameter', 'medium_index', 'medium_material', 'index', 'material', pre=True)
    def validate_array(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    @validator('diameter', pre=True)
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

        self.binding = CppCylinderSet(**binding_kwargs)

# -
