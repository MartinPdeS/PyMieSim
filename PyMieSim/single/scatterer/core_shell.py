#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyvista
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, RefractiveIndex
from PyMieSim.single.source.base import BaseSource

from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_scatterer import CORESHELL
from PyMieSim.single.source.base import BaseSource


class CoreShell(CORESHELL, BaseScatterer):
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

    def __init__(self,
        core_diameter: Length,
        shell_thickness: Length,
        core_property: RefractiveIndex | BaseMaterial,
        shell_property: RefractiveIndex | BaseMaterial,
        medium_property: RefractiveIndex | BaseMaterial,
        source: BaseSource,
    ) -> None:
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
        core_diameter = Length.check(core_diameter)
        shell_thickness = Length.check(shell_thickness)
        source = BaseSource.check(source)

        core_index, self.core_material = self._assign_index_or_material(
            core_property
        )
        shell_index, self.shell_material = self._assign_index_or_material(
            shell_property
        )
        medium_index, self.medium_material = self._assign_index_or_material(
            medium_property
        )

        super().__init__(
            core_diameter=core_diameter,
            shell_thickness=shell_thickness,
            core_refractive_index=core_index,
            shell_refractive_index=shell_index,
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
        sphape = pyvista.Sphere(
            center=(0.0, 0.0, 0.0), radius=0.1, theta_resolution=100, phi_resolution=100
        )

        # Add the cone mesh to the scene
        scene.add_mesh(sphape, color=color, opacity=opacity)
