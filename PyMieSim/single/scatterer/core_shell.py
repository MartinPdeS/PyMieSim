#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pyvista
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, Volume, RefractiveIndex, ureg
from pydantic.dataclasses import dataclass

from PyMieSim.utils import config_dict
from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_scatterer import CORESHELL
from PyMieSim.single.source.base import BaseSource


@dataclass(config=config_dict, kw_only=True)
class CoreShell(CORESHELL, BaseScatterer):
    """
    Class representing a core/shell spherical scatterer.

    Parameters
    ----------
    core_diameter : Length
        Diameter of the core of the scatterer.
    shell_thickness : Length
        Thickness of the shell surrounding the core.
    core_property : RefractiveIndex | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer's core. Only one can be provided.
    shell_property : RefractiveIndex | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer's shell. Only one can be provided.
    medium_property : ureg.Quantity or BaseMaterial
        Refractive index or material of the surrounding medium.
    source : BaseSource
        Source object associated with the scatterer.
    """

    core_diameter: Length
    shell_thickness: Length
    core_property: RefractiveIndex | BaseMaterial
    shell_property: RefractiveIndex | BaseMaterial
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
        Initialize the CoreShell scatterer with its core and shell properties.

        """
        self.core_index, self.core_material = self._assign_index_or_material(
            self.core_property
        )
        self.shell_index, self.shell_material = self._assign_index_or_material(
            self.shell_property
        )
        self.medium_index, self.medium_material = self._assign_index_or_material(
            self.medium_property
        )

        super().__init__(
            core_diameter=self.core_diameter.to(ureg.meter).magnitude,
            shell_thickness=self.shell_thickness.to(ureg.meter).magnitude,
            core_refractive_index=self.core_index.to(ureg.RIU).magnitude,
            shell_refractive_index=self.shell_index.to(ureg.RIU).magnitude,
            medium_refractive_index=self.medium_index.to(ureg.RIU).magnitude,
            source=self.source,
        )

    @property
    def radius(self) -> Length:
        """Return the outer radius of the scatterer."""
        return self.core_diameter / 2 + self.shell_thickness

    @property
    def volume(self) -> Volume:
        """Return the volume of the scatterer."""
        vol = (4 / 3) * numpy.pi * (self.radius**3)
        return vol.to(ureg.meter**3)

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
        sphape = pyvista.Sphere(
            center=(0.0, 0.0, 0.0), radius=0.1, theta_resolution=100, phi_resolution=100
        )

        # Add the cone mesh to the scene
        scene.add_mesh(sphape, color=color, opacity=opacity)
