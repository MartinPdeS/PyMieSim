#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.units import ureg
from PyMieSim.single.mesh import FibonacciMesh


@pytest.fixture
def fibonacci_mesh():
    """
    Fixture to create a FibonacciMesh instance with predefined parameters.

    Returns
    -------
    FibonacciMesh
        Instance of FibonacciMesh with fixed test parameters.
    """
    return FibonacciMesh(
        max_angle=np.pi / 4 * ureg.radian,  # Maximum angle (45 degrees in radians)
        min_angle=0 * ureg.radian,  # Minimum angle (0 degrees in radians)
        sampling=100,  # Number of sample points
        phi_offset=30 * ureg.degree,  # Offset in the azimuthal angle (in degrees)
        gamma_offset=45 * ureg.degree,  # Offset in the polar angle (in degrees)
        rotation=0 * ureg.degree,  # Rotation angle (in degrees)
    )


def test_initialization(fibonacci_mesh):
    """
    Test that the FibonacciMesh object is initialized with the correct parameters.
    """
    value = np.pi / 4 * ureg.radian
    assert fibonacci_mesh.max_angle == value
    assert fibonacci_mesh.sampling == 100
    assert np.isclose(fibonacci_mesh.phi_offset.to('degree').magnitude, 30)
    assert np.isclose(fibonacci_mesh.gamma_offset.to('degree').magnitude, 45)
    assert np.isclose(fibonacci_mesh.rotation.to('degree').magnitude, 0)


def test_projection_HV_vector(fibonacci_mesh):
    """
    Test the projection of vectors into parallel and perpendicular components.
    """
    parallel_projection, perpendicular_projection = (
        fibonacci_mesh.projection_HV_vector()
    )

    # Ensure the shape of the projections is as expected
    assert parallel_projection.shape == (2, fibonacci_mesh.sampling, 3)
    assert perpendicular_projection.shape == (2, fibonacci_mesh.sampling, 3)


def test_projection_HV_scalar(fibonacci_mesh):
    """
    Test the scalar projection into parallel and perpendicular components.
    """
    parallel_projection, perpendicular_projection = (
        fibonacci_mesh.projection_HV_scalar()
    )

    # Ensure the shape of the scalar projections is as expected
    assert parallel_projection.shape == (2, fibonacci_mesh.sampling)
    assert perpendicular_projection.shape == (2, fibonacci_mesh.sampling)


def test_get_cartesian_coordinates(fibonacci_mesh):
    """
    Test that Cartesian coordinates are correctly computed and returned.
    """
    coordinates = np.vstack(
        [
            fibonacci_mesh.cartesian.x,
            fibonacci_mesh.cartesian.y,
            fibonacci_mesh.cartesian.z,
        ]
    )

    # Ensure the shape of the Cartesian coordinates is as expected
    assert coordinates.shape == (3, fibonacci_mesh.sampling)

def test_rotate_around_axis(fibonacci_mesh):
    """
    Test that the mesh is correctly rotated around its axis.
    """
    initial_coordinates = np.vstack(
        [
            fibonacci_mesh.cartesian.x,
            fibonacci_mesh.cartesian.y,
            fibonacci_mesh.cartesian.z,
        ]
    )

    # Rotate the mesh around its axis by 45 degrees
    fibonacci_mesh.rotate_around_axis(45)

    rotate_coordinates = np.vstack(
        [
            fibonacci_mesh.cartesian.x,
            fibonacci_mesh.cartesian.y,
            fibonacci_mesh.cartesian.z,
        ]
    )

    # Ensure the coordinates have changed after the rotation
    assert not np.allclose(initial_coordinates, rotate_coordinates, atol=1e-6)


if __name__ == "__main__":
    pytest.main(["-W error", "-s", __file__])
