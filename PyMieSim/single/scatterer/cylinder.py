#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, RefractiveIndex

from PyMieSim.single.source.base import BaseSource
from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_single import CYLINDER


class Cylinder(CYLINDER, BaseScatterer):
    property_names = [
        "size_parameter",
        "radius",
        "cross_section",
        "g",
        "Qsca",
        "Qext",
        "Qabs",
        "Csca",
        "Cext",
        "Cabs",
    ]

    def __init__(
        self,
        diameter: Length,
        property: RefractiveIndex | BaseMaterial,
        medium_property: RefractiveIndex | BaseMaterial,
        source: BaseSource,
    ) -> None:
        """
        Represents a cylindrical scatterer used for scattering simulations in optical systems.

        Parameters
        ----------
        diameter : units.Quantity
            Diameter of the cylindrical scatterer, given in meters.
        property : units.Quantity | BaseMaterial
            Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer. Only one of these should be provided at a time to specify the core characteristic.
        medium_property : RefractiveIndex or BaseMaterial
            Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the surrounding medium. Only one can be provided.
        source : Gaussian | PlaneWave
            The source object associated with the scatterer.
        """
        diameter = Length.check(diameter)
        source = BaseSource.check(source)

        index, self.material = self._assign_index_or_material(
            wavelength=source.wavelength,
            property=property
        )
        medium_index, self.medium_material = self._assign_index_or_material(
            wavelength=source.wavelength,
            property=medium_property
        )

        super().__init__(
            diameter=diameter,
            refractive_index=index,
            medium_refractive_index=medium_index,
            source=source,
        )
