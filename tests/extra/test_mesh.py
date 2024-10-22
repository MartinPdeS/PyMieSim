#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.mesh import FibonacciMesh


@pytest.fixture
def fibonacci_mesh():
    """
    Fixture to create a FibonacciMesh instance with predefined parameters.

    Returns:
        FibonacciMesh: Instance of FibonacciMesh with fixed test parameters.
    """
    return FibonacciMesh(
        max_angle=np.pi / 4,   # Maximum angle (45 degrees in radians)
        sampling=1000,         # Number of sample points
        phi_offset=30,         # Offset in the azimuthal angle (in degrees)
        gamma_offset=45,       # Offset in the polar angle (in degrees)
        rotation_angle=0       # Rotation angle (in degrees)
    )


def test_initialization(fibonacci_mesh):
    """
    Test that the FibonacciMesh object is initialized with the correct parameters.
    """
    assert fibonacci_mesh.max_angle == np.pi / 4
    assert fibonacci_mesh.sampling == 1000
    assert fibonacci_mesh.phi_offset == 30
    assert fibonacci_mesh.gamma_offset == 45
    assert fibonacci_mesh.rotation_angle == 0


def test_get_phi(fibonacci_mesh):
    """
    Test that the azimuthal angles (phi) are returned in both radians and degrees correctly.
    """
    phi_radians = fibonacci_mesh.get_phi(unit='radian')
    phi_degrees = fibonacci_mesh.get_phi(unit='degree')

    # Ensure phi values are consistent between radians and degrees
    assert np.allclose(phi_radians, np.deg2rad(phi_degrees), atol=1e-6)


def test_get_theta(fibonacci_mesh):
    """
    Test that the polar angles (theta) are returned in both radians and degrees correctly.
    """
    theta_radians = fibonacci_mesh.get_theta(unit='radian')
    theta_degrees = fibonacci_mesh.get_theta(unit='degree')

    # Ensure theta values are consistent between radians and degrees
    assert np.allclose(theta_radians, np.deg2rad(theta_degrees), atol=1e-6)


def test_projection_HV_vector(fibonacci_mesh):
    """
    Test the projection of vectors into parallel and perpendicular components.
    """
    parallel_projection, perpendicular_projection = fibonacci_mesh.projection_HV_vector()

    # Ensure the shape of the projections is as expected
    assert parallel_projection.shape == (2, fibonacci_mesh.sampling, 3)
    assert perpendicular_projection.shape == (2, fibonacci_mesh.sampling, 3)


def test_projection_HV_scalar(fibonacci_mesh):
    """
    Test the scalar projection into parallel and perpendicular components.
    """
    parallel_projection, perpendicular_projection = fibonacci_mesh.projection_HV_scalar()

    # Ensure the shape of the scalar projections is as expected
    assert parallel_projection.shape == (2, fibonacci_mesh.sampling)
    assert perpendicular_projection.shape == (2, fibonacci_mesh.sampling)


def test_get_cartesian_coordinates(fibonacci_mesh):
    """
    Test that Cartesian coordinates are correctly computed and returned.
    """
    coordinates = fibonacci_mesh.get_cartesian_coordinates()

    # Ensure the shape of the Cartesian coordinates is as expected
    assert coordinates.shape == (3, fibonacci_mesh.sampling)


def test_get_axis_vector(fibonacci_mesh):
    """
    Test that the axis vector is correctly computed and has unit length.
    """
    axis_vector = fibonacci_mesh.get_axis_vector()

    # Ensure the axis vector has the correct shape and is normalized to length 1
    assert axis_vector.shape == (3,)
    assert np.isclose(np.linalg.norm(axis_vector), 1.0, atol=1e-6)


def test_rotate_around_axis(fibonacci_mesh):
    """
    Test that the mesh is correctly rotated around its axis.
    """
    initial_coordinates = fibonacci_mesh.get_cartesian_coordinates()

    # Rotate the mesh around its axis by 45 degrees
    fibonacci_mesh.rotate_around_axis(45)

    rotated_coordinates = fibonacci_mesh.get_cartesian_coordinates()

    # Ensure the coordinates have changed after the rotation
    assert not np.allclose(initial_coordinates, rotated_coordinates, atol=1e-6)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
