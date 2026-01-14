#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyvista
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, RefractiveIndex

from PyMieSim.single.source.base import BaseSource
from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_scatterer import CYLINDER


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

        index, self.material = self._assign_index_or_material(property)
        medium_index, self.medium_material = self._assign_index_or_material(medium_property)

        super().__init__(
            diameter=diameter,
            refractive_index=index,
            medium_refractive_index=medium_index,
            source=source,
        )

    def _add_to_3d_ax(
        self, scene: pyvista.Plotter, color: str = "black", opacity: float = 1.0
    ) -> None:
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
            resolution=100,  # Number of sides for the cylinder
        )

        # Add the cone mesh to the scene
        scene.add_mesh(shape, color=color, opacity=opacity)
