#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyvista
from PyOptik.material.base_class import BaseMaterial

from PyMieSim import units
from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_scatterer import CORESHELL
from PyMieSim.single.source.base import BaseSource


class CoreShell(CORESHELL, BaseScatterer):
    """
    Class representing a core/shell spherical scatterer.

    Parameters
    ----------
    core_diameter : units.Quantity
        Diameter of the core of the scatterer.
    shell_thickness : units.Quantity
        Thickness of the shell surrounding the core.
    core_property : units.Quantity | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer's core. Only one can be provided.
    shell_property : units.Quantity | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer's shell. Only one can be provided.

    """

    core_diameter: units.Quantity
    shell_thickness: units.Quantity
    core_property: units.Quantity | BaseMaterial
    shell_property: units.Quantity | BaseMaterial
    medium_property: units.Quantity | BaseMaterial
    source: BaseSource

    property_names = [
        "size_parameter", "cross_section", "g",
        "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qpr",
        "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cpr"
    ]

    def __init__(self,
        core_diameter: units.Quantity,
        shell_thickness: units.Quantity,
        core_property: units.Quantity | BaseMaterial,
        shell_property: units.Quantity | BaseMaterial,
        medium_property: units.Quantity | BaseMaterial,
        source: BaseSource):
        """
        Initialize the CoreShell scatterer with its core and shell properties.

        Parameters
        ----------
        core_diameter : units.Quantity
            Diameter of the core in meters.
        shell_thickness : units.Quantity
            Thickness of the shell in meters.
        core_property : units.Quantity or BaseMaterial
            Refractive index or material of the core.
        shell_property : units.Quantity or BaseMaterial
            Refractive index or material of the shell.
        medium_property : units.Quantity or BaseMaterial
            Refractive index or material of the surrounding medium.
        source : BaseSource
            Source object associated with the scatterer.
        """
        self.core_diameter = self._validate_units(core_diameter, dimension='distance', units=units.meter)
        self.shell_thickness = self._validate_units(shell_thickness, dimension='distance', units=units.meter)
        self.core_property = self._validate_property(core_property)
        self.shell_property = self._validate_property(shell_property)
        self.medium_property = self._validate_property(medium_property)

        self.source = source
        self.core_index, self.core_material = self._assign_index_or_material(self.core_property)
        self.shell_index, self.shell_material = self._assign_index_or_material(self.shell_property)
        self.medium_index, self.medium_material = self._assign_index_or_material(self.medium_property)

        super().__init__(
            core_diameter=core_diameter.to(units.meter).magnitude,
            shell_thickness=shell_thickness.to(units.meter).magnitude,
            core_refractive_index=self.core_index.to(units.RIU).magnitude,
            shell_refractive_index=self.shell_index.to(units.RIU).magnitude,
            medium_refractive_index=self.medium_index.to(units.RIU).magnitude,
            source=self.source
        )

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
        sphape = pyvista.Sphere(
            center=(0.0, 0.0, 0.0),
            radius=0.1,
            theta_resolution=100,
            phi_resolution=100
        )

        # Add the cone mesh to the scene
        scene.add_mesh(sphape, color=color, opacity=opacity)