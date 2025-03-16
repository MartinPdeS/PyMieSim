#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union
import numpy
from pydantic.dataclasses import dataclass
from dataclasses import field
from PyMieSim.polarization import BasePolarization
from PyMieSim.special_functions import NA_to_angle
from PyMieSim.binary.interface_source import BindedGaussian # type: ignore
import pyvista
from PyMieSim.units import Quantity
from PyMieSim.single.source.base import BaseSource, config_dict


@dataclass(config=config_dict)
class Gaussian(BaseSource):
    """
    Represents a Gaussian light source for light scattering simulations, characterized by its optical power and numerical aperture.

    Parameters
    ----------
    wavelength: Quantity
        Wavelength of the light field in meters.
    polarization : BasePolarization | Quantity
        Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
    optical_power: Quantity
        Optical power of the source in Watts.
    NA: Quantity
        Numerical aperture of the source.
    """
    wavelength: Quantity
    polarization: Union[BasePolarization, Quantity]
    optical_power: Quantity
    NA: Quantity

    amplitude: Quantity = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.polarization, BasePolarization):
            self.polarization = BasePolarization(self.polarization)

        self.wavenumber = 2 * numpy.pi / self.wavelength
        self.waist = self.wavelength / (self.NA * numpy.pi)
        self.peak_intensity = 2 * self.optical_power / (numpy.pi * self.waist**2)

        self.binding = BindedGaussian(
            wavelength=self.wavelength.to_base_units().magnitude,
            NA=self.NA.to_base_units().magnitude,
            optical_power=self.optical_power.to_base_units().magnitude,
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

        # Add a translucent sphere to the scene
        sphere = pyvista.Sphere(radius=1)
        scene.add_mesh(sphere, opacity=0.3)

        # Display the scene
        scene.show()

    def _add_to_3d_ax(self, scene: pyvista.Plotter, color: str = 'red', opacity: float = 0.8) -> None:
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
        # Calculate the maximum angle from the numerical aperture (NA)

        max_angle = numpy.rad2deg(NA_to_angle(NA=self.NA.magnitude))

        # Create the cone mesh
        cone_mesh = pyvista.Cone(
            center=(0.0, 0.0, -0.5),
            direction=(0.0, 0.0, 1.0),
            height=0.9,
            resolution=100,
            angle=max_angle
        )

        # Add the cone mesh to the scene
        scene.add_mesh(cone_mesh, color=color, opacity=0.3)
# -
