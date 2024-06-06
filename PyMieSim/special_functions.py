#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


def get_rotation_matrix(axis: numpy.ndarray, theta: float) -> numpy.ndarray:
    """
    Calculate the rotation matrix for a counterclockwise rotation around a given axis.

    Parameters:
        - axis: numpy.ndarray, the axis of rotation.
        - theta: float, the angle of rotation in degrees.

    Returns:
        - numpy.ndarray, the rotation matrix.
    """
    theta_rad = numpy.deg2rad(theta)
    axis = axis / numpy.linalg.norm(axis)
    a, b, c = axis * numpy.sin(theta_rad / 2.0)
    aa, bb, cc, dd = numpy.cos(theta_rad / 2.0), a * a, b * b, c * c
    bc, ad, ac, ab, bd, cd = b * c, a * numpy.cos(theta_rad / 2.0), a * c, a * b, b * numpy.cos(theta_rad / 2.0), c * numpy.cos(theta_rad / 2.0)

    return numpy.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def NA_to_angle(NA: float) -> float:
    """
    Convert numerical aperture (NA) to angle in radians.

    Parameters:
        - NA: float, the numerical aperture.

    Returns:
        - float, the angle in radians.
    """
    return numpy.arcsin(NA) if NA <= 1.0 else numpy.arcsin(NA - 1) + numpy.pi / 2


def angle_space_to_direct(angle_space: numpy.ndarray, k: float) -> numpy.ndarray:
    """
    Convert angle space to direct space.

    Parameters:
        - angle_space: numpy.ndarray, the angle space array.
        - k: float, the wave number.

    Returns:
    - numpy.ndarray, the direct space array.
    """
    rad_space = numpy.deg2rad(angle_space)
    fourier_space = numpy.sin(rad_space) * k / (2 * numpy.pi)
    fourier_unit = numpy.abs(numpy.diff(fourier_space)[0])
    return numpy.fft.fftshift(numpy.fft.fftfreq(angle_space.size, d=fourier_unit))


def direct_space_to_angle(direct_space: numpy.ndarray, k: float) -> numpy.ndarray:
    """
    Convert direct space to angle space.

    Parameters:
        - direct_space: numpy.ndarray, the direct space array.
        - k: float, the wave number.

    Returns:
    - numpy.ndarray, the angle space array in degrees.
    """
    direct_unit = numpy.abs(numpy.diff(direct_space)[0])
    fourier_space = numpy.fft.fftshift(numpy.fft.fftfreq(direct_space.size, d=direct_unit))
    angle_space = numpy.arcsin(2 * numpy.pi * fourier_space / k)
    if numpy.isnan(angle_space).any():
        raise ValueError("Magnification too large.")
    return numpy.rad2deg(angle_space)


def cartesian_to_spherical(x: numpy.ndarray, y: numpy.ndarray, z: numpy.ndarray) -> tuple:
    """
    Convert Cartesian coordinates to spherical coordinates.

    Parameters:
        - x: numpy.ndarray, the x coordinates.
        - y: numpy.ndarray, the y coordinates.
        - z: numpy.ndarray, the z coordinates.

    Returns:
        - tuple of numpy.ndarray, (r, phi, theta) the spherical coordinates.
    """
    r = numpy.sqrt(x**2 + y**2 + z**2)
    phi = numpy.arcsin(z / r)
    theta = numpy.arctan2(y, x)
    return r, phi, theta


def spherical_to_cartesian(phi: numpy.ndarray, theta: numpy.ndarray, r: numpy.ndarray = None) -> tuple:
    """
    Convert spherical coordinates to Cartesian coordinates.

    Parameters:
        - phi: numpy.ndarray, the phi angles.
        - theta: numpy.ndarray, the theta angles.
        - r: numpy.ndarray, the radial distances; defaults to unit radius if None.

    Returns:
        - tuple of numpy.ndarray, (x, y, z) the Cartesian coordinates.
    """
    if r is None:
        r = numpy.ones_like(phi)

    x = r * numpy.cos(phi) * numpy.cos(theta)
    y = r * numpy.cos(phi) * numpy.sin(theta)
    z = r * numpy.sin(phi)
    return x, y, z


def rotate_on_y(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    """
    Rotate spherical coordinates around the Y-axis.

    Parameters:
        - phi: numpy.ndarray, the phi angles in radians.
        - theta: numpy.ndarray, the theta angles in radians.
        - angle: float, the rotation angle in radians.

    Returns:
        - tuple of numpy.ndarray, the rotated spherical coordinates (r, phi, theta).
    """
    # Convert to Cartesian coordinates for rotation
    x, y, z = spherical_to_cartesian(phi, theta)
    # Apply rotation around the Y-axis
    xp = x * numpy.cos(angle) + z * numpy.sin(angle)
    zp = -x * numpy.sin(angle) + z * numpy.cos(angle)
    # Convert back to spherical coordinates
    return cartesian_to_spherical(xp, y, zp)


def rotate_on_z(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    """
    Rotate spherical coordinates around the Z-axis.

    Parameters:
        - phi: numpy.ndarray, the phi angles in radians.
        - theta: numpy.ndarray, the theta angles in radians.
        - angle: float, the rotation angle in radians.

    Returns:
        - tuple of numpy.ndarray, the rotated spherical coordinates (r, phi, theta).
    """
    # Convert to Cartesian for rotation
    x, y, z = spherical_to_cartesian(phi, theta)
    # Apply rotation around the Z-axis
    xp = x * numpy.cos(angle) - y * numpy.sin(angle)
    yp = x * numpy.sin(angle) + y * numpy.cos(angle)
    # Convert back to spherical coordinates
    return cartesian_to_spherical(xp, yp, z)


def rotate_on_x(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    """
    Rotate spherical coordinates around the X-axis.

    Parameters:
        - phi: numpy.ndarray, the phi angles in radians.
        - theta: numpy.ndarray, the theta angles in radians.
        - angle: float, the rotation angle in radians.

    Returns:
        - tuple of numpy.ndarray, the rotated spherical coordinates (r, phi, theta).
    """
    # Convert to Cartesian for rotation
    x, y, z = spherical_to_cartesian(phi, theta)
    # Apply rotation around the X-axis
    yp = y * numpy.cos(angle) - z * numpy.sin(angle)
    zp = y * numpy.sin(angle) + z * numpy.cos(angle)
    # Convert back to spherical coordinates
    return cartesian_to_spherical(x, yp, zp)


# -
