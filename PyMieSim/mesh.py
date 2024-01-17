#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from dataclasses import dataclass

from PyMieSim.binary.Fibonacci import FIBONACCI
from MPSPlots.render3D import SceneList as SceneList3D


@dataclass
class FibonacciMesh():
    """
    Class wich represent an angular mesh. The distribution of points inside
    the mesh is similar to a Fibonacci sphere where each point cover an
    equivalent solid angle.
    """
    max_angle: float = 1.5
    """ Angle in radian defined by the numerical aperture of the imaging system. """
    sampling: int = 1000
    """Number of point distrubuted inside the Solid angle defined by the numerical aperture. """
    phi_offset: float = 0.
    """ Angle offset in the parallel direction of the polarization of incindent light in degree. """
    gamma_offset: float = 0.
    """ Angle offset in the perpendicular direction of the polarization of incindent light in degree. """
    rotation_angle: float = 0
    """ Rotation of the mesh around it's principal axis. """

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
        self.generate_ledvedev_mesh()

    @property
    def plan(self):
        if self._plan is None:
            self.cpp_binding.Computeplan()
            self._plan = numpy.asarray([self.cpp_binding.Xplan, self.cpp_binding.Yplan, self.cpp_binding.Zplan])

        return self._plan

    @property
    def perpendicular_vector(self):
        if self._perp is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.perpVector

    @property
    def parallel_vector(self):
        if self._para is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.paraVector

    @property
    def horizontal_to_perpendicular(self):
        if self._horizontal_to_perpendicular_vector is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.horizontal_to_perpendicular_vector

    @property
    def horizontal_to_parallel(self):
        if self._horizontal_to_parallel_vector is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.horizontal_to_parallel_vector

    @property
    def vertical_to_perpendicular(self):
        if self._vertical_to_perpendicular_vector is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.perpVector

    @property
    def vertical_to_parallel(self):
        if self._vertical_to_parallel_vector is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.paraVector

    @property
    def parallel_plan(self):
        if self._parallel_vector_in_plan is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.paraVectorZplanar

    @property
    def perpendicular_plan(self):
        if self._perpendicular_vector_in_plan is None:
            self.cpp_binding.compute_vector_field()

        return self.cpp_binding.perpVectorZplanar

    def get_phi(self, unit: str = 'radian'):
        if unit.lower() in ['rad', 'radian']:
            return self.cpp_binding.phi

        elif unit.lower() in ['deg', 'degree']:
            return numpy.rad2deg(self.cpp_binding.phi)

    def get_theta(self, unit: str = 'radian'):
        if unit.lower() in ['rad', 'radian']:
            return self.cpp_binding.theta

        elif unit.lower() in ['deg', 'degree']:
            return numpy.rad2deg(self.cpp_binding.theta)

    @property
    def X(self):
        return self.cpp_binding.x

    @property
    def Y(self):
        return self.cpp_binding.y

    @property
    def Z(self):
        return self.cpp_binding.z

    def initialize_properties(self):
        self.cartesian_coordinates = numpy.c_[
            self.cpp_binding.x,
            self.cpp_binding.y,
            self.cpp_binding.z
        ].T

        self.d_omega_radian = self.cpp_binding.d_omega
        self.d_omega_degree = self.cpp_binding.d_omega * (180 / numpy.pi)**2

        self.omega_radian = self.cpp_binding.omega
        self.omega_degree = self.cpp_binding.omega * (180 / numpy.pi)**2

    def projection_HV_vector(self) -> tuple:
        parallel_projection = [
            self.projection_on_base_vector(
                vector=self.vertical_to_parallel_plan,
                base_vector=X
            ) for X in [self.vertical_vector, self.horizontal_vector]
        ]

        perpendicular_projection = [
            self.projection_on_base_vector(
                vector=self.vertical_to_perpendicular_plan,
                base_vector=X
            ) for X in [self.vertical_vector, self.horizontal_vector]
        ]

        return numpy.array(parallel_projection), numpy.array(perpendicular_projection)

    def projection_HV_scalar(self) -> tuple:
        parallel_projection = [
            self.projection_on_base_scalar(
                vector=self.vertical_to_perpendicular_in_z_plan,
                base_vector=X
            ) for X in [self.vertical_vector, self.horizontal_vector]
        ]

        perpendicular_projection = [
            self.projection_on_base_scalar(
                vector=self.vertical_to_perpendicular_in_z_plan,
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

    def plot(self) -> None:
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

    def _add_mesh_to_ax_(self, ax) -> None:
        axis_vector = self.get_axis_vector()

        coordinates = self.get_cartesian_coordinates()

        ax.add_unstructured_mesh(coordinates=coordinates)

        ax.add_unstructured_mesh(coordinates=axis_vector * 1.4)

        # ax.add_unit_sphere()
        ax.add_unit_axis()
        ax.add_unit_theta_vector(radius=1.1)
        ax.add_unit_phi_vector(radius=1.1)

    def rotate_around_axis(self, rotation_angle: float) -> None:
        """
        Rotate the mesh around its principal axis.

        :param      rotation_angle:  The rotation angle [degree]
        :type       rotation_angle:  float

        :returns:   No returns
        :rtype:     None
        """
        self.cpp_binding.rotate_around_axis(rotation_angle)

    def generate_ledvedev_mesh(self) -> None:
        self.cpp_binding = FIBONACCI(
            sampling=self.sampling,
            max_angle=self.max_angle,
            phi_offset=numpy.deg2rad(self.phi_offset),
            gamma_offset=numpy.deg2rad(self.gamma_offset),
            rotation_angle=self.rotation_angle
        )

        self.initialize_properties()


# -
