#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from typing import Optional
from pydantic.dataclasses import dataclass
from TypedUnit import Angle, ureg

from PyMieSim.binary.interface_detector import FIBONACCIMESH
from PyMieSim.utils import config_dict


@dataclass(config=config_dict, kw_only=True)
class FibonacciMesh(FIBONACCIMESH):
    """
    Represents an angular mesh using a Fibonacci sphere distribution, where each point
    covers an equivalent solid angle. This type of mesh is useful for simulating angular
    distributions in light scattering experiments.

    The mesh is defined by a numerical aperture and a number of sampling points, with
    additional control over the orientation and rotation of the mesh.

    Parameters
    ----------
    max_angle : Angle
        The maximum angle in radians, typically defined by the numerical aperture of the imaging system.
    sampling : int
        The number of points distributed across the mesh. Higher values result in finer resolution.
    phi_offset : Angle
        Angle offset in the azimuthal (parallel to incident light polarization) direction in degrees.
    gamma_offset : Angle
        Angle offset in the polar (perpendicular to incident light polarization) direction in degrees.
    rotation_angle : Optional[Angle], default=0.0
        Rotation of the entire mesh around its principal axis, in degrees.

    Attributes
    ----------
    cartesian_coordinates : numpy.ndarray
        Cartesian coordinates of the mesh points, calculated during initialization.
    d_omega_radian : float
        Differential solid angle element in radians for each point on the mesh.
    d_omega_degree : float
        Differential solid angle element in degrees for each point on the mesh.
    omega_radian : float
        Total solid angle covered by the mesh in radians.
    omega_degree : float
        Total solid angle covered by the mesh in degrees.

    Methods
    -------
    projection_HV_vector() -> tuple:
        Returns the parallel and perpendicular projections of the mesh on the horizontal and vertical base vectors as vectors.
    projection_HV_scalar() -> tuple:
        Returns the parallel and perpendicular projections of the mesh on the horizontal and vertical base vectors as scalar values.
    projection_on_base_scalar(vector, base_vector) -> numpy.ndarray:
        Computes the scalar projection of a given vector onto a base vector.
    projection_on_base_vector(vector, base_vector) -> numpy.ndarray:
        Computes the vector projection of a given vector onto a base vector.
    get_axis_vector() -> numpy.ndarray:
        Returns the axis vector corresponding to the phi and gamma offsets.
    rotate_around_axis(rotation_angle) -> None:
        Rotates the mesh around its principal axis by a specified angle.

    """
    max_angle: Angle
    sampling: int
    phi_offset: Angle
    gamma_offset: Angle
    min_angle: Angle = 0 * ureg.degree
    rotation_angle: Optional[Angle] = 0.0 * ureg.degree

    # ---------------------- Initialization Methods ----------------------
    def __post_init__(self):
        """
        Initializes the FibonacciMesh with the specified parameters.
        """
        super().__init__(
            sampling=self.sampling,
            max_angle=self.max_angle,
            min_angle=self.min_angle,
            rotation_angle=self.rotation,
            phi_offset=self.phi_offset.to(ureg.radian).magnitude,
            gamma_offset=self.gamma_offset.to(ureg.radian).magnitude,
        )

        self.compute_vector_field()

    # ---------------------- Property Methods ----------------------
    @property
    def theta(self) -> Angle:
        """
        Returns the polar angles (theta) of the mesh in radians.

        Returns
        -------
        units.Quantity
            The polar angles in radians.
        """
        return self.spherical.theta * ureg.radian

    @property
    def phi(self) -> Angle:
        """
        Returns the azimuthal angles (phi) of the mesh in radians.

        Returns
        -------
        units.Quantity
            The azimuthal angles in radians.
        """
        return self.spherical.phi * ureg.radian

    @property
    def d_omega(self) -> Angle:
        """
        Returns the solid angle (d_omega) of the mesh in steradians.

        Returns
        -------
        units.Quantity
            The solid angle in steradians.
        """
        return self.spherical._cpp_d_omega * ureg.steradian

    @property
    def omega(self) -> Angle:
        """
        Returns the solid angle (omega) of the mesh in steradians.

        Returns
        -------
        units.Quantity
            The solid angle in steradians.
        """
        return self.spherical._cpp_omega * ureg.steradian

    @property
    def rotation(self) -> Angle:
        """
        Returns the rotation angle of the mesh around its principal axis in radians.

        Returns
        -------
        units.Quantity
            The rotation angle in radians.
        """
        return self._cpp_rotation * ureg.radian

    # ---------------------- Additional Methods ----------------------
    def get_axis_vector(self) -> numpy.ndarray:
        """
        Returns the axis vector corresponding to the phi and gamma offsets.

        Returns
        -------
        numpy.ndarray
            The axis vector normalized to unit length.
        """
        axis_vector = numpy.array([self.cartesian.x[0], self.cartesian.y[0], self.cartesian.z[0]])

        norm = numpy.sqrt(numpy.square(axis_vector).sum())

        return axis_vector / norm

    def projection_HV_vector(self) -> tuple:
        """
        Projects the vertical and horizontal base vectors onto the parallel and perpendicular components
        of the mesh, returning these projections as vectors.

        The projections help visualize the relationship between the mesh's angular distribution and
        the incident light polarization directions.

        Returns
        -------
        tuple of numpy.ndarray
            A tuple where the first element is the projection of the vertical and horizontal base vectors
            onto the parallel vector field, and the second element is the projection onto the perpendicular
            vector field.
            Each of these projections is a 2D array of shape (2, 3), where:
                - The first row corresponds to the projection onto the vertical base vector.
                - The second row corresponds to the projection onto the horizontal base vector.
        """
        parallel_projection = [
            self.projection_on_base_vector(
                vector=self.parallel,
                base_vector=X
            ) for X in [self._cpp_vertical_base, self._cpp_horizontal_base]
        ]

        perpendicular_projection = [
            self.projection_on_base_vector(
                vector=self.perpendicular,
                base_vector=X
            ) for X in [self._cpp_vertical_base, self._cpp_horizontal_base]
        ]

        return numpy.array(parallel_projection), numpy.array(perpendicular_projection)

    def projection_HV_scalar(self) -> tuple:
        """
        Projects the vertical and horizontal base vectors onto the parallel and perpendicular components
        of the mesh, returning these projections as scalar values (dot products).

        This method calculates the scalar projection, which represents the magnitude of the component
        of the vector field along the base vectors.

        Returns
        -------
        tuple of numpy.ndarray
            A tuple where the first element is the scalar projection of the vertical and horizontal base vectors
            onto the parallel vector field, and the second element is the scalar projection onto the perpendicular
            vector field.
            Each of these projections is a 1D array of length 2, where:
                - The first element corresponds to the projection onto the vertical base vector.
                - The second element corresponds to the projection onto the horizontal base vector.

        """
        parallel_projection = [
            self.projection_on_base_scalar(
                vector=self.parallel,
                base_vector=X
            ) for X in [self._cpp_vertical_base, self._cpp_horizontal_base]
        ]

        perpendicular_projection = [
            self.projection_on_base_scalar(
                vector=self.perpendicular,
                base_vector=X
            ) for X in [self._cpp_vertical_base, self._cpp_horizontal_base]
        ]

        return numpy.array(parallel_projection), numpy.array(perpendicular_projection)

    def projection_on_base_scalar(self, vector: numpy.ndarray, base_vector: numpy.ndarray) -> numpy.ndarray:
        r"""
        Computes the scalar projection (dot product) of a given vector onto a base vector.

        This method calculates the magnitude of the projection of `vector` along the direction
        of `base_vector`. The scalar projection is useful for understanding how much of
        `vector` is aligned with `base_vector`.

        Parameters
        ----------
        vector : numpy.ndarray
            The vector to be projected. It should be a 1D array representing the vector field component
            (either parallel or perpendicular) of the mesh.
        base_vector : numpy.ndarray
            The base vector onto which the projection is made. It should be a 1D array, typically
            representing the horizontal or vertical base direction.

        Returns
        -------
        numpy.ndarray
            The scalar projection of `vector` onto `base_vector`, which is a single numeric value.
            This value represents the magnitude of the component of `vector` that is aligned with
            `base_vector`.

        Notes
        -----
        The scalar projection is calculated as the dot product of `vector` and `base_vector`:

        .. math::
            \text{scalar projection} = \vec{v} \cdot \vec{b}

        where :math:`\vec{v}` is the input vector and :math:`\vec{b}` is the base vector.
        """
        return vector._cpp_get_scalar_product(base_vector)

    def projection_on_base_vector(self, vector: numpy.ndarray, base_vector: numpy.ndarray) -> numpy.ndarray:
        r"""
        Computes the vector projection of a given vector onto a base vector.

        This method calculates the projection of `vector` along the direction of `base_vector`
        and returns a new vector that represents the component of `vector` aligned with
        `base_vector`. This is useful for decomposing vectors in terms of directional components.

        Parameters
        ----------
        vector : numpy.ndarray
            The vector to be projected. It should be a 1D array representing the vector field component
            (either parallel or perpendicular) of the mesh.
        base_vector : numpy.ndarray
            The base vector onto which the projection is made. It should be a 1D array, typically
            representing the horizontal or vertical base direction.

        Returns
        -------
        numpy.ndarray
            The vector projection of `vector` onto `base_vector`, which is a vector in the direction
            of `base_vector`. The resulting vector represents the part of `vector` that lies along
            `base_vector`.

        Notes
        -----
        The vector projection is calculated as:

        .. math::
            \text{vector projection} = \left( \frac{\vec{v} \cdot \vec{b}}{\vec{b} \cdot \vec{b}} \right) \vec{b}

        where :math:`\vec{v}` is the input vector and :math:`\vec{b}` is the base vector. This formula gives the projection
        of `vector` along the direction of `base_vector`.

        The method first computes the scalar projection of `vector` onto `base_vector` using the dot product, and
        then scales `base_vector` by this scalar to produce the projected vector.
        """
        projection = self.projection_on_base_scalar(vector, base_vector)
        base_projection = numpy.outer(projection, base_vector.data)

        return base_projection
