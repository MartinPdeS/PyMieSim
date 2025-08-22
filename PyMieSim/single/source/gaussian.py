#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pyvista
from TypedUnit import Length, ElectricField, ureg, Angle, Power, Dimensionless
from pydantic.dataclasses import dataclass

from PyMieSim.utils import config_dict
from PyMieSim.polarization import BasePolarization
from PyMieSim.binary.interface_source import GAUSSIAN
from PyMieSim.single.source.base import BaseSource


@dataclass(config=config_dict, kw_only=True)
class Gaussian(GAUSSIAN, BaseSource):
    """
    Represents a Gaussian light source for optical simulations.

    wavelength: Length
        Wavelength of the light field in meters.
    polarization : Angle | BasePolarization
        Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
    optical_power: Power
        Optical power of the source in Watts.
    NA: Dimensionless
        Numerical aperture of the source.
    """
    wavelength: Length
    polarization: Angle | BasePolarization
    optical_power: Power
    NA: Dimensionless

    amplitude: ElectricField = None

    def __post_init__(self) -> None:
        """
        Initializes a Gaussian light source for light scattering simulations, characterized by its optical power and numerical aperture.

        """
        self.wavenumber = 2 * numpy.pi / self.wavelength
        self.waist = self.wavelength / (self.NA * numpy.pi)
        self.peak_intensity = 2 * self.optical_power / (numpy.pi * self.waist**2)
        self.polarization = self._validate_source_polarization(self.polarization)


        super().__init__(
            wavelength=self.wavelength.to(ureg.meter).magnitude,
            jones_vector=self.polarization.element[0],
            optical_power=self.optical_power.to(ureg.watt).magnitude,
            NA=self.NA.to(ureg.AU).magnitude,
        )

    def numerical_aperture_to_angle(self, numerical_aperture: float) -> float:
        """
        Convert numerical aperture (NA) to angle in radians.

        Parameters
        ----------
        NA : float
            The numerical aperture.

        Returns
        -------
        float
            The angle in radians.
        """
        return numpy.arcsin(numerical_aperture) if numerical_aperture <= 1.0 else numpy.arcsin(numerical_aperture - 1) + numpy.pi / 2

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

        numerical_aperture = self.NA.magnitude

        angle = self.numerical_aperture_to_angle(numerical_aperture)

        max_angle = numpy.rad2deg(angle)

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
