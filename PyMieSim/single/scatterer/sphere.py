#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pyvista
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, Volume, RefractiveIndex, ureg
from pydantic.dataclasses import dataclass

from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.utils import config_dict
from PyMieSim.binary.interface_scatterer import SPHERE
from PyMieSim.single.source.base import BaseSource


@dataclass(config=config_dict, kw_only=True)
class Sphere(SPHERE, BaseScatterer):
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

    diameter: Length
    property: RefractiveIndex | BaseMaterial
    medium_property: RefractiveIndex | BaseMaterial
    source: BaseSource

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

    def __post_init__(self):
        """
        Initialize the Sphere scatterer with its diameter and material properties.
        """
        self.index, self.material = self._assign_index_or_material(self.property)
        self.medium_index, self.medium_material = self._assign_index_or_material(
            self.medium_property
        )

        super().__init__(
            diameter=self.diameter.to(ureg.meter).magnitude,
            refractive_index=self.index.to(ureg.RIU).magnitude,
            medium_refractive_index=self.medium_index.to(ureg.RIU).magnitude,
            source=self.source,
        )

    @property
    def radius(self) -> Length:
        """Return the radius of the sphere."""
        return self.diameter / 2

    @property
    def volume(self) -> Volume:
        """Return the volume of the sphere."""
        return (4 / 3) * numpy.pi * (self.radius**3)

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
        shape = pyvista.Sphere(
            center=(0.0, 0.0, 0.0), radius=0.1, theta_resolution=100, phi_resolution=100
        )

        # Add the cone mesh to the scene
        scene.add_mesh(shape, color=color, opacity=opacity)
