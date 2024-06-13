import pytest
import numpy as np
from PyMieSim.mesh import FibonacciMesh


@pytest.fixture
def fibonacci_mesh():
    return FibonacciMesh(
        max_angle=np.pi / 4,  # 45 degrees in radians
        sampling=1000,
        phi_offset=30,  # degrees
        gamma_offset=45,  # degrees
        rotation_angle=0  # degrees
    )


def test_initialization(fibonacci_mesh):
    assert fibonacci_mesh.max_angle == np.pi / 4
    assert fibonacci_mesh.sampling == 1000
    assert fibonacci_mesh.phi_offset == 30
    assert fibonacci_mesh.gamma_offset == 45
    assert fibonacci_mesh.rotation_angle == 0


def test_get_phi(fibonacci_mesh):
    phi_radians = fibonacci_mesh.get_phi(unit='radian')
    phi_degrees = fibonacci_mesh.get_phi(unit='degree')
    assert np.allclose(phi_radians, np.deg2rad(phi_degrees), atol=1e-6)


def test_get_theta(fibonacci_mesh):
    theta_radians = fibonacci_mesh.get_theta(unit='radian')
    theta_degrees = fibonacci_mesh.get_theta(unit='degree')
    assert np.allclose(theta_radians, np.deg2rad(theta_degrees), atol=1e-6)


def test_projection_HV_vector(fibonacci_mesh):
    parallel_projection, perpendicular_projection = fibonacci_mesh.projection_HV_vector()
    assert parallel_projection.shape == (2, fibonacci_mesh.sampling, 3)
    assert perpendicular_projection.shape == (2, fibonacci_mesh.sampling, 3)


def test_projection_HV_scalar(fibonacci_mesh):
    parallel_projection, perpendicular_projection = fibonacci_mesh.projection_HV_scalar()
    assert parallel_projection.shape == (2, fibonacci_mesh.sampling)
    assert perpendicular_projection.shape == (2, fibonacci_mesh.sampling)


def test_get_cartesian_coordinates(fibonacci_mesh):
    coordinates = fibonacci_mesh.get_cartesian_coordinates()
    assert coordinates.shape == (3, fibonacci_mesh.sampling)


def test_get_axis_vector(fibonacci_mesh):
    axis_vector = fibonacci_mesh.get_axis_vector()
    assert axis_vector.shape == (3,)
    assert np.isclose(np.linalg.norm(axis_vector), 1.0, atol=1e-6)


def test_rotate_around_axis(fibonacci_mesh):
    initial_coordinates = fibonacci_mesh.get_cartesian_coordinates()
    fibonacci_mesh.rotate_around_axis(45)  # Rotate by 45 degrees
    rotated_coordinates = fibonacci_mesh.get_cartesian_coordinates()
    assert not np.allclose(initial_coordinates, rotated_coordinates, atol=1e-6)


if __name__ == "__main__":
    pytest.main()
