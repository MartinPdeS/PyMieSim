#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Optional, NoReturn
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict
import numpy

from PyMieSim.binary.Fibonacci import FibonacciMesh as CPPFibonacciMesh
from MPSPlots.render3D import SceneList as SceneList3D

config_dict = ConfigDict(
    kw_only=True,
    extra='forbid',
    slots=True
)


@dataclass(config=config_dict)
class FibonacciMesh:
    """
    Represents an angular mesh using a Fibonacci sphere distribution where each point covers
    an equivalent solid angle. This setup is designed to simulate angular distributions
    in light scattering experiments.
    """
    max_angle: float  # Angle in radians defined by the numerical aperture of the imaging system.
    sampling: int  # Number of points distributed inside the solid angle.
    phi_offset: float  # Angle offset in the parallel direction of incident light polarization in degrees.
    gamma_offset: float  # Angle offset in the perpendicular direction of incident light polarization in degrees.
    rotation_angle: Optional[float] = 0.  # Rotation of the mesh around its principal axis in degrees.

    def __post_init__(self):
        self.structured = False

        self._para = None
        self._perp = None

        self._parallel_vector_in_plan = None
        self._perpendicular_vector_in_plan = None

        self._vertical_to_perpendicular_vector = None
        self._horizontal_to_perpendicular_vector = None
        self._vertical_to_parallel_vector = None
        self._horizontal_to_parallel_vector = None

        self._phi = None
        self._theta = None

        self._plan = None

        self.vertical_vector = numpy.array([1, 0, 0])
        self.horizontal_vector = numpy.array([0, 1, 0])

        self.binding = CPPFibonacciMesh(
            sampling=self.sampling,
            max_angle=self.max_angle,
            phi_offset=numpy.deg2rad(self.phi_offset),
            gamma_offset=numpy.deg2rad(self.gamma_offset),
            rotation_angle=self.rotation_angle
        )

        self.initialize_properties()

        self.binding.compute_vector_field()

    @property
    def perpendicular_vector(self):
        return self.binding.perpendicular_vector

    @property
    def parallel_vector(self):
        return self.binding.parallel_vector

    @property
    def horizontal_to_perpendicular(self):
        return self.binding.horizontal_to_perpendicular_vector

    @property
    def horizontal_to_parallel(self):
        return self.binding.horizontal_to_parallel_vector

    @property
    def vertical_to_perpendicular(self):
        return self.binding.perpendicular_vector

    @property
    def vertical_to_parallel(self):
        return self.binding.parallel_vector

    def get_phi(self, unit: str = 'radian'):
        if unit.lower() in ['rad', 'radian']:
            return self.binding.phi

        elif unit.lower() in ['deg', 'degree']:
            return numpy.rad2deg(self.binding.phi)

    def get_theta(self, unit: str = 'radian'):
        if unit.lower() in ['rad', 'radian']:
            return self.binding.theta

        elif unit.lower() in ['deg', 'degree']:
            return numpy.rad2deg(self.binding.theta)

    @property
    def X(self) -> numpy.ndarray:
        return self.binding.x

    @property
    def Y(self) -> numpy.ndarray:
        return self.binding.y

    @property
    def Z(self) -> numpy.ndarray:
        return self.binding.z

    def initialize_properties(self):
        """
        Initializes additional mesh properties based on the Fibonacci distribution
        generated by the binding.
        """
        self.cartesian_coordinates = numpy.column_stack((self.binding.x, self.binding.y, self.binding.z))

        self.d_omega_radian = self.binding.d_omega
        self.d_omega_degree = self.binding.d_omega * (180 / numpy.pi)**2

        self.omega_radian = self.binding.omega
        self.omega_degree = self.binding.omega * (180 / numpy.pi)**2

    def projection_HV_vector(self) -> tuple:
        parallel_projection = [
            self.projection_on_base_vector(
                vector=self.vertical_to_parallel,
                base_vector=X
            ) for X in [self.vertical_vector, self.horizontal_vector]
        ]

        perpendicular_projection = [
            self.projection_on_base_vector(
                vector=self.vertical_to_perpendicular,
                base_vector=X
            ) for X in [self.vertical_vector, self.horizontal_vector]
        ]

        return numpy.array(parallel_projection), numpy.array(perpendicular_projection)

    def projection_HV_scalar(self) -> tuple:
        parallel_projection = [
            self.projection_on_base_scalar(
                vector=self.vertical_to_parallel,
                base_vector=X
            ) for X in [self.vertical_vector, self.horizontal_vector]
        ]

        perpendicular_projection = [
            self.projection_on_base_scalar(
                vector=self.vertical_to_perpendicular,
                base_vector=X
            ) for X in [self.vertical_vector, self.horizontal_vector]
        ]

        return numpy.array(parallel_projection), numpy.array(perpendicular_projection)

    def projection_on_base_scalar(self, vector: numpy.ndarray, base_vector: numpy.ndarray) -> numpy.ndarray:
        return vector.dot(base_vector)

    def projection_on_base_vector(self, vector: numpy.ndarray, base_vector: numpy.ndarray) -> numpy.ndarray:
        projection = self.projection_on_base_scalar(vector, base_vector)

        base_projection = numpy.outer(projection, base_vector)

        return base_projection

    def plot(self) -> SceneList3D:
        """
        Plots the Fibonacci mesh using 3D rendering.

        Returns:
            SceneList3D: The 3D scene list object containing the mesh plot.
        """
        figure = SceneList3D()
        ax = figure.append_ax()
        self._add_mesh_to_ax_(ax)

        return figure

    def get_cartesian_coordinates(self) -> numpy.ndarray:
        return numpy.c_[self.X, self.Y, self.Z].T

    def get_axis_vector(self) -> numpy.ndarray:
        """
        Returns the axis vector that correspond ot the phi and gamma offset

        :returns:   The axis vector.
        :rtype:     numpy.ndarray
        """
        x, y, z = self.get_cartesian_coordinates()
        axis_vector = numpy.array([x[0], y[0], z[0]])

        norm = numpy.sqrt(numpy.square(axis_vector).sum())

        return axis_vector / norm

    def _add_mesh_to_ax_(self, ax) -> NoReturn:
        axis_vector = self.get_axis_vector()

        coordinates = self.get_cartesian_coordinates()

        ax.add_unstructured_mesh(coordinates=coordinates)

        ax.add_unstructured_mesh(coordinates=axis_vector * 1.4)

        # ax.add_unit_sphere()
        ax.add_unit_axis()
        ax.add_unit_theta_vector(radius=1.1)
        ax.add_unit_phi_vector(radius=1.1)

    def rotate_around_axis(self, rotation_angle: float) -> NoReturn:
        """
        Rotate the mesh around its principal axis.

        :param      rotation_angle:  The rotation angle [degree]
        :type       rotation_angle:  float

        :returns:   No returns
        :rtype:     None
        """
        self.binding.rotate_around_axis(rotation_angle)


# -
