#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import List
from pydantic.dataclasses import dataclass
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, RefractiveIndex, AnyUnit

from PyMieSim.binary.interface_sets import CppSphereSet
from PyMieSim.experiment.scatterer.base import BaseScatterer
from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import Sequential
from PyMieSim.utils import config_dict


@dataclass(config=config_dict, kw_only=True)
class Sphere(BaseScatterer, Sequential):
    """
    A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

    This class provides specific implementations for setting up and binding spherical scatterers
    with their properties to a simulation environment. It extends the `BaseScatterer` class by
    adding spherical-specific attributes and methods for handling simulation setups.

    Parameters
    ----------
    source : PyMieSim.experiment.source.base.BaseSource
        Light source configuration for the simulation.
    diameter : Length
        Diameter(s) of the spherical scatterers in meters.
    property : List[BaseMaterial] | List[RefractiveIndex]
        Refractive index or indices of the spherical scatterers themselves.
    medium_property : List, optional
        BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
    """

    source: BaseSource
    diameter: Length
    property: List[BaseMaterial] | List[RefractiveIndex]
    medium_property: List[BaseMaterial] | List[RefractiveIndex]

    available_measure_list = [
        "Qsca",
        "Qext",
        "Qabs",
        "Qratio",
        "Qforward",
        "Qback",
        "Qpr",
        "Csca",
        "Cext",
        "Cabs",
        "Cratio",
        "Cforward",
        "Cback",
        "Cpr",
        "a1",
        "a2",
        "a3",
        "b1",
        "b2",
        "b3",
        "g",
        "coupling",
    ]

    attributes = ["diameter", "property", "medium_property"]

    def _generate_binding(self):
        """
        Constructs the keyword arguments necessary for the C++ binding interface, specifically tailored for spherical scatterers.
        This includes processing material indices and organizing them into a structured dictionary for simulation interaction.

        This method automatically configures the `binding_kwargs` attribute with appropriately formatted values.
        """
        self.mapping = {}

        medium_properties = self._add_properties(
            name="medium", properties=self.medium_property
        )

        scatterer_properties = self._add_properties(
            name="scatterer", properties=self.property
        )

        self.binding_kwargs = dict(
            diameter=self.diameter,
            is_sequential=self.is_sequential,
            medium_properties=medium_properties,
            scatterer_properties=scatterer_properties,
        )

        self.set = CppSphereSet(
            **{
                k: v.to_base_units().magnitude if isinstance(v, AnyUnit) else v
                for k, v in self.binding_kwargs.items()
            }
        )
