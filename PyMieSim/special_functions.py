#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


def NA_to_angle(NA: float) -> float:
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
    return numpy.arcsin(NA) if NA <= 1.0 else numpy.arcsin(NA - 1) + numpy.pi / 2


def cartesian_to_spherical(x: numpy.ndarray, y: numpy.ndarray, z: numpy.ndarray) -> tuple:
    """
    Convert Cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    x : numpy.ndarray
        The x coordinates.
    y : numpy.ndarray
        The y coordinates.
    z : numpy.ndarray
        The z coordinates.

    Returns
    -------
    numpy.ndarray
        The spherical coordinates (r, phi, theta).
    """
    r = numpy.sqrt(x**2 + y**2 + z**2)
    phi = numpy.arcsin(z / r)
    theta = numpy.arctan2(y, x)
    return r, phi, theta


def spherical_to_cartesian(phi: numpy.ndarray, theta: numpy.ndarray, r: numpy.ndarray = None) -> tuple:
    """
    Convert spherical coordinates to Cartesian coordinates.

    Parameters
    ----------
    phi : numpy.ndarray
        The phi angles.
    theta : numpy.ndarray
        The theta angles.
    r : numpy.ndarray
        The radial distances; defaults to unit radius if None.

    Returns
    -------
    numpy.ndarray
        The Cartesian coordinates (x, y, z).
    """
    if r is None:
        r = numpy.ones_like(phi)

    x = r * numpy.cos(phi) * numpy.cos(theta)
    y = r * numpy.cos(phi) * numpy.sin(theta)
    z = r * numpy.sin(phi)
    return x, y, z


def rotate_on_x(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    """
    Rotate spherical coordinates around the X-axis.

    Parameters
    ----------
    phi : numpy.ndarray
        The phi angles.
    theta : numpy.ndarray
        The theta angles.
    r : numpy.ndarray
        The radial distances; defaults to unit radius if None.

    Returns
    -------
    numpy.ndarray
        The Cartesian coordinates (x, y, z).
    """
    # Convert to Cartesian for rotation
    x, y, z = spherical_to_cartesian(phi, theta)
    # Apply rotation around the X-axis
    yp = y * numpy.cos(angle) - z * numpy.sin(angle)
    zp = y * numpy.sin(angle) + z * numpy.cos(angle)
    # Convert back to spherical coordinates
    return cartesian_to_spherical(x, yp, zp)
