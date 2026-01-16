#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, RefractiveIndex

from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_single import SPHERE
from PyMieSim.single.source.base import BaseSource


class Sphere(SPHERE, BaseScatterer):
    """
    Class representing a homogeneous spherical scatterer.
    """
    property_names = [
        "size_parameter",
        "radius",
        "volume",
        "cross_section",
        "g",
        "Qsca",
        "Qext",
        "Qabs",
        "Qback",
        "Qratio",
        "Qpr",
        "Csca",
        "Cext",
        "Cabs",
        "Cback",
        "Cratio",
        "Cpr",
    ]

    def __init__(
        self,
        diameter: Length,
        property: RefractiveIndex | BaseMaterial,
        medium_property: RefractiveIndex | BaseMaterial,
        source: BaseSource,
    ) -> None:
        """
        Class representing a homogeneous spherical scatterer.

        Parameters
        ----------
        diameter : Length
            Diameter of the cylindrical scatterer, given in meters.
        property : RefractiveIndex or BaseMaterial
            Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer. Only one can be provided.
        medium_property : RefractiveIndex or BaseMaterial
            Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the surrounding medium. Only one can be provided.
        source : BaseSource
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
        print("sphere: ", medium_index)
        super().__init__(
            diameter=diameter,
            refractive_index=index,
            medium_refractive_index=medium_index,
            source=source,
        )