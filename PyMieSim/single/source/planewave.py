#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union
import numpy
from pydantic.dataclasses import dataclass
from PyMieSim.polarization import BasePolarization
from PyMieSim.binary.interface_source import BindedPlanewave # type: ignore
import pyvista
from PyMieSim.units import Quantity
from PyMieSim.single.source.base import BaseSource, config_dict


@dataclass(config=config_dict)
class PlaneWave(BaseSource):
    """
    Represents a plane wave light source for light scattering simulations.

    Inherits from LightSource and specifies amplitude directly.

    Parameters
    ----------
    wavelength : Quantity
        Wavelength of the light field in meters.
    polarization : BasePolarization | Quantity
        Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
    amplitude : Quantity
        Amplitude of the electric field.

    """
    wavelength: Quantity
    amplitude: Quantity
    polarization: Union[BasePolarization, Quantity]

    def __post_init__(self) -> None:
        if not isinstance(self.polarization, BasePolarization):
            self.polarization = BasePolarization(self.polarization)

        self.wavenumber = 2 * numpy.pi / self.wavelength

        self.binding = BindedPlanewave(
            wavelength=self.wavelength.to_base_units().magnitude,
            amplitude=self.amplitude.to_base_units().magnitude,
            jones_vector=self.polarization.element[0]
        )

    def plot(self, color: str = 'red', opacity: float = 0.8, show_axis_label: bool = False) -> None:
        """
        Plots the 3D structure of the Gaussian source.

        This method creates a 3D plot of the Gaussian source, adds the structure to the plot,
        and optionally displays axis labels.

        Parameters
        ----------
        color : str
            The color of the structure in the plot. Default is 'red'.
        opacity : float
            The opacity of the structure. Default is 0.8.
        show_axis_label : bool
            If True, axis labels will be shown. Default is False.

        """
        # Create a 3D plotting scene
        scene = pyvista.Plotter()

        # Add the structure to the scene
        self._add_to_3d_ax(scene=scene, color=color, opacity=opacity)

        # Add axes at the origin, optionally showing axis labels
        scene.add_axes_at_origin(labels_off=not show_axis_label)

        # Display the scene
        scene.show()

    def _add_to_3d_ax(self, scene: pyvista.Plotter, color: str = 'red', opacity: float = 0.8) -> None:
        """
        Adds a 3D cone representation to the given PyVista plotting scene.

        The cylinder represents the acceptance angle determined by the numerical aperture (NA) of the system.
        The cylinder is positioned at the origin and points downward along the z-axis.

        Parameters
        ----------
        scene : pyvista.Plotter
            The 3D plotting scene to which the cone will be added.
        color : str
            The color of the cone mesh. Default is 'red'.
        opacity : float
            The opacity of the cone mesh. Default is 0.8.

        """
        # Define the cylinder parameters
        cylinder_mesh = pyvista.Cylinder(
            center=(0.0, 0.0, 0.5),
            direction=(0.0, 0.0, -1.0),
            radius=0.2,
            height=1.0,
            resolution=100
        )

        # Add the cone mesh to the scene
        scene.add_mesh(cylinder_mesh, color='blue', opacity=0.3)
