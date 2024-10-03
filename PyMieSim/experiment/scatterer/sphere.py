#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from pydantic import field_validator
from typing import List

import numpy
from PyMieSim.binary.SetsInterface import CppSphereSet
from PyOptik.base_class import BaseMaterial
from PyMieSim.units import Quantity, meter, RIU
from PyMieSim.experiment.scatterer.base_class import BaseScatterer, config_dict

@dataclass(config=config_dict, kw_only=True)
class Sphere(BaseScatterer):
    """
    A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

    This class provides specific implementations for setting up and binding spherical scatterers
    with their properties to a simulation environment. It extends the `BaseScatterer` class by
    adding spherical-specific attributes and methods for handling simulation setups.

    Attributes:
        source (PyMieSim.experiment.source.base.BaseSource): Light source configuration for the simulation.
        diameter (Quantity): Diameter(s) of the spherical scatterers in meters.
        property (List[BaseMaterial] | List[Quantity]): Refractive index or indices of the spherical scatterers themselves.
        medium_property (List, optional): BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
    """
    diameter: Quantity
    property: List[BaseMaterial] | List[Quantity]

    available_measure_list = [
        'Qsca', 'Qext', 'Qabs', 'Qratio', 'Qforward', 'Qback', 'Qpr',
        'Csca', 'Cext', 'Cabs', 'Cratio', 'Cforward', 'Cback', 'Cpr', 'a1',
        'a2', 'a3', 'b1', 'b2', 'b3', 'g', 'coupling',
    ]

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

        self._assign_index_or_material(element='', property=self.property)
        self._assign_index_or_material(element='medium_', property=self.medium_property)

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppSphereSet(**binding_kwargs)