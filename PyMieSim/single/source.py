#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union, NoReturn
import numpy
from pydantic.dataclasses import dataclass
from dataclasses import field

from PyMieSim.polarization import UnitPolarizationAngle
from PyMieSim.special_functions import NA_to_angle
from MPSPlots.render3D import SceneList as SceneList3D
from PyMieSim.binary.SourceInterface import BindedGaussian, BindedPlanewave, BindedBaseSource  # noqa: F401
import pyvista

config_dict = dict(
    kw_only=True,
    slots=True,
    extra='forbid'
)


@dataclass(config=config_dict)
class PlaneWave():
    """
    Represents a plane wave light source for light scattering simulations.

    Inherits from LightSource and specifies amplitude directly.

    Attributes:
        wavelength (float): Wavelength of the light field in meters.
        polarization (Union[UnitPolarizationAngle, float]): Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        amplitude (float): Amplitude of the electric field.
    """
    wavelength: float
    polarization: Union[UnitPolarizationAngle, float]
    amplitude: float = 1

    def __post_init__(self) -> NoReturn:
        if not isinstance(self.polarization, UnitPolarizationAngle):
            self.polarization = UnitPolarizationAngle(self.polarization)

        self.wavenumber = 2 * numpy.pi / self.wavelength

        self.binding = BindedPlanewave(
            wavelength=self.wavelength,
            amplitude=self.amplitude,
            jones_vector=self.polarization.jones_vector
        )

    def plot(self, color: str = 'red', opacity: float = 0.8, show_axis_label: bool = False) -> NoReturn:
        """
        Plots the 3D structure of the Gaussian source.

        This method creates a 3D plot of the Gaussian source, adds the structure to the plot,
        and optionally displays axis labels.

        Args:
            color (str): The color of the structure in the plot. Default is 'red'.
            opacity (float): The opacity of the structure. Default is 0.8.
            show_axis_label (bool): If True, axis labels will be shown. Default is False.

        Returns:
            NoReturn: This method does not return a value. It displays the 3D plot.
        """
        # Create a 3D plotting scene
        scene = pyvista.Plotter()

        # Add the structure to the scene
        self._add_to_3d_ax(scene=scene, color=color, opacity=opacity)

        # Add axes at the origin, optionally showing axis labels
        scene.add_axes_at_origin(labels_off=not show_axis_label)

        # Display the scene
        scene.show()


    def _add_to_3d_ax(self, scene: pyvista.Plotter, color: str = 'red', opacity: float = 0.8) -> NoReturn:
        """
        Adds a 3D cone representation to the given PyVista plotting scene.

        The cylinder represents the acceptance angle determined by the numerical aperture (NA) of the system.
        The cylinder is positioned at the origin and points downward along the z-axis.

        Args:
            scene (pyvista.Plotter): The 3D plotting scene to which the cone will be added.
            color (str): The color of the cone mesh. Default is 'red'.
            opacity (float): The opacity of the cone mesh. Default is 0.8.

        Returns:
            NoReturn: This method does not return a value. It adds the cone mesh to the provided scene.
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


@dataclass(config=config_dict)
class Gaussian():
    """
    Represents a Gaussian light source for light scattering simulations, characterized by its optical power and numerical aperture.

    Attributes:
        wavelength (float): Wavelength of the light field in meters.
        polarization_value (float): Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        amplitude (float): Amplitude of the electric field.
        optical_power (float): Optical power of the source in Watts.
        NA (float): Numerical aperture of the source.
    """
    wavelength: float
    polarization: Union[UnitPolarizationAngle, float]
    optical_power: float
    NA: float
    amplitude: float = field(init=False, repr=False)

    def __post_init__(self) -> NoReturn:
        if not isinstance(self.polarization, UnitPolarizationAngle):
            self.polarization = UnitPolarizationAngle(self.polarization)

        self.wavenumber = 2 * numpy.pi / self.wavelength

        self.binding = BindedGaussian(
            wavelength=self.wavelength,
            NA=self.NA,
            optical_power=self.optical_power,
            jones_vector=self.polarization.jones_vector
        )

    def plot(self, color: str = 'red', opacity: float = 0.8, show_axis_label: bool = False) -> NoReturn:
        """
        Plots the 3D structure of the Gaussian source.

        This method creates a 3D plot of the Gaussian source, adds the structure to the plot,
        and optionally displays axis labels.

        Args:
            color (str): The color of the structure in the plot. Default is 'red'.
            opacity (float): The opacity of the structure. Default is 0.8.
            show_axis_label (bool): If True, axis labels will be shown. Default is False.

        Returns:
            NoReturn: This method does not return a value. It displays the 3D plot.
        """
        # Create a 3D plotting scene
        scene = pyvista.Plotter()

        # Add the structure to the scene
        self._add_to_3d_ax(scene=scene, color=color, opacity=opacity)

        # Add axes at the origin, optionally showing axis labels
        scene.add_axes_at_origin(labels_off=not show_axis_label)

        # Display the scene
        scene.show()


    def _add_to_3d_ax(self, scene: pyvista.Plotter, color: str = 'red', opacity: float = 0.8) -> NoReturn:
        """
        Adds a 3D cone representation to the given PyVista plotting scene.

        The cone represents the acceptance angle determined by the numerical aperture (NA) of the system.
        The cone is positioned at the origin and points downward along the z-axis.

        Args:
            scene (pyvista.Plotter): The 3D plotting scene to which the cone will be added.
            color (str): The color of the cone mesh. Default is 'red'.
            opacity (float): The opacity of the cone mesh. Default is 0.8.

        Returns:
            NoReturn: This method does not return a value. It adds the cone mesh to the provided scene.
        """
        # Calculate the maximum angle from the numerical aperture (NA)
        max_angle = numpy.rad2deg(NA_to_angle(NA=self.NA))

        # Create the cone mesh
        cone_mesh = pyvista.Cone(
            center=(0.0, 0.0, 0.45),
            direction=(0.0, 0.0, -1.0),
            height=0.9,
            resolution=100,
            angle=max_angle
        )

        # Add the cone mesh to the scene
        scene.add_mesh(cone_mesh, color='blue', opacity=0.3)
# -
