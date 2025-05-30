#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyvista
from PyOptik.material.base_class import BaseMaterial

from PyMieSim import units
from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_scatterer import CYLINDER
from PyMieSim.single.source.base import BaseSource

class Cylinder(CYLINDER, BaseScatterer):
    """
    Represents a cylindrical scatterer used for scattering simulations in optical systems.

    Parameters
    ----------
    diameter : units.Quantity
        Diameter of the cylindrical scatterer, given in meters.
    property : units.Quantity | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer. Only one of these should be provided at a time to specify the core characteristic.

    """
    diameter: units.Quantity
    property: units.Quantity | BaseMaterial
    medium_property: units.Quantity | BaseMaterial
    source: BaseSource

    property_names = [
        "size_parameter", "cross_section", "g",
        "Qsca", "Qext", "Qabs",
        "Csca", "Cext", "Cabs"
    ]

    def __init__(self, diameter: units.Quantity, property: units.Quantity | BaseMaterial, medium_property: units.Quantity | BaseMaterial, source: BaseSource):
        """
        Initialize the Cylinder scatterer with its diameter and material properties.

        Parameters
        ----------
        diameter : units.Quantity
            Diameter of the sphere in meters.
        property : units.Quantity or BaseMaterial
            Refractive index or material of the sphere.
        medium_property : units.Quantity
            Refractive index of the surrounding medium.
        source : Source
            Source object associated with the scatterer.
        """
        self.diameter = self._validate_units(diameter, dimension='distance', units=units.meter)
        self.property = self._validate_property(property)
        self.medium_property = self._validate_property(medium_property)
        self.source = source

        self.index, self.material = self._assign_index_or_material(self.property)
        self.medium_index, self.medium_material = self._assign_index_or_material(self.medium_property)

        super().__init__(
            diameter=diameter.to(units.meter).magnitude,
            refractive_index=self.index.to(units.RIU).magnitude,
            medium_refractive_index=self.medium_index.to(units.RIU).magnitude,
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

    def _add_to_3d_ax(self, scene: pyvista.Plotter, color: str = 'black', opacity: float = 1.0) -> None:
        """
        Adds a 3D cone representation to the given PyVista plotting scene.

        The cone represents the acceptance angle determined by the numerical aperture (NA) of the system.
        The cone is positioned at the origin and points downward along the z-axis.

        Parameters
        ----------
        scene : pyvista.Plotter
            The 3D plotting scene to which the cone will be added.
        color : str
            The color of the cone mesh. Default is 'red'.
        opacity : float
            The opacity of the cone mesh. Default is 0.8.

        """
        # Create the cone mesh
        shape = pyvista.Cylinder(
            center=(0.0, 0.0, 0.0),
            radius=0.1,
            height=2.0,  # Height of the cylinder
            direction=(0, -1, 0),  # Pointing downwards along the z-axis
            resolution=100  # Number of sides for the cylinder
        )

        # Add the cone mesh to the scene
        scene.add_mesh(shape, color=color, opacity=opacity)