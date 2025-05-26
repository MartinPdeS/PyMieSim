#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import Quantity
from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_scatterer import CYLINDER
from PyMieSim.single.source.base import BaseSource

class Cylinder(CYLINDER, BaseScatterer):
    """
    Represents a cylindrical scatterer used for scattering simulations in optical systems.

    Parameters
    ----------
    diameter : Quantity
        Diameter of the cylindrical scatterer, given in meters.
    property : Quantity | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer. Only one of these should be provided at a time to specify the core characteristic.

    """
    diameter: Quantity
    property: Quantity | BaseMaterial
    medium_property: Quantity | BaseMaterial
    source: BaseSource

    property_names = [
        "size_parameter", "cross_section", "g",
        "Qsca", "Qext", "Qabs",
        "Csca", "Cext", "Cabs"
    ]

    def __init__(self, diameter: Quantity, property: Quantity | BaseMaterial, medium_property: Quantity | BaseMaterial, source: BaseSource):
        """
        Initialize the Cylinder scatterer with its diameter and material properties.

        Parameters
        ----------
        diameter : Quantity
            Diameter of the sphere in meters.
        property : Quantity or BaseMaterial
            Refractive index or material of the sphere.
        medium_property : Optional[Quantity]
            Refractive index of the surrounding medium.
        source : Optional[Source]
            Source object associated with the scatterer.
        """
        self.diameter = self._validate_length(diameter)
        self.property = self._validate_property(property)
        self.medium_property = self._validate_property(medium_property)
        self.source = source

        self.index, self.material = self._assign_index_or_material(self.property)
        self.medium_index, self.medium_material = self._assign_index_or_material(self.medium_property)

        super().__init__(
            diameter=diameter.to_base_units().magnitude,
            refractive_index=self.index.to_base_units().magnitude,
            medium_refractive_index=self.medium_index.to_base_units().magnitude,
            source=self.source
        )

    @property
    def Cback(self) -> None:
        raise NotImplementedError

    @property
    def Qback(self) -> None:
        raise NotImplementedError

    @property
    def Cratio(self) -> None:
        raise NotImplementedError

    @property
    def Qratio(self) -> None:
        raise NotImplementedError
