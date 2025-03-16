#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List
from PyMieSim.units import Quantity
from pydantic.dataclasses import dataclass
from PyMieSim.binary.interface_sets import CppCylinderSet
from PyOptik.material.base_class import BaseMaterial
from PyMieSim.experiment.scatterer.base import BaseScatterer
from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import config_dict, Sequential


@dataclass(config=config_dict, kw_only=True)
class Cylinder(BaseScatterer, Sequential):
    """
    Represents a cylindrical scatterer configuration for PyMieSim simulations.

    Parameters
    ----------
    source : PyMieSim.experiment.source.base.BaseSource
        Light source configuration for the simulation.
    diameter : Quantity
        Diameter(s) of the cylinder in meters.
    property : List[BaseMaterial] | List[Quantity]
        Refractive index or indices of the spherical scatterers themselves.
    medium_property : List[BaseMaterial] | List[Quantity]
        BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
    """
    source: BaseSource
    diameter: Quantity
    property: List[BaseMaterial] | List[Quantity]
    medium_property: List[BaseMaterial] | List[Quantity]

    available_measure_list = [
        'Qsca', 'Qext', 'Qabs', 'Csca', 'Cext', 'Cabs', 'a11',
        'a21', 'a12', 'a22', 'a13', 'a23', 'b11', 'b21', 'b12',
        'b22', 'b13', 'b23', 'coupling',
    ]

    attributes = ['diameter', 'property', 'medium_property']

    def _generate_binding(self) -> None:
        """
        Constructs the keyword arguments necessary for the C++ binding interface, specifically tailored for spherical scatterers.
        This includes processing material indices and organizing them into a structured dictionary for simulation interaction.

        This method automatically configures the `binding_kwargs` attribute with appropriately formatted values.
        """
        self.mapping = {}

        self.binding_kwargs = dict(diameter=self.diameter, is_sequential=self.is_sequential)

        self._add_properties(name='medium', properties=self.medium_property)

        self._add_properties(name='scatterer', properties=self.property)

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppCylinderSet(**binding_kwargs)
