#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy

# Configuration dictionary for the Pydantic dataclass

def cartesian_to_spherical(
    x: numpy.ndarray, y: numpy.ndarray, z: numpy.ndarray
) -> tuple:
    r"""Convert Cartesian coordinates to spherical coordinates.

    PyMieSim uses the elevation angle ``phi`` and the azimuth ``theta``:

    .. math::

        r = \sqrt{x^2 + y^2 + z^2}, \qquad
        \phi = \arcsin\left(\frac{z}{r}\right), \qquad
        \theta = \operatorname{atan2}(y, x).

    The returned angles are in radians and follow the same convention as
    :func:`spherical_to_cartesian`. The inputs are broadcast by NumPy.

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
    tuple of numpy.ndarray
        The spherical coordinates ``(r, phi, theta)``. ``r`` has the input
        length unit, while ``phi`` and ``theta`` are NumPy arrays containing
        radians.
    """
    r = numpy.sqrt(x**2 + y**2 + z**2)
    phi = numpy.arcsin(z / r)
    theta = numpy.arctan2(y, x)
    return r, phi, theta


def spherical_to_cartesian(
    phi: numpy.ndarray, theta: numpy.ndarray, r: numpy.ndarray = None
) -> tuple:
    r"""Convert spherical coordinates to Cartesian coordinates.

    ``phi`` is the elevation from the :math:`x-y` plane and ``theta`` is the
    azimuth in that plane. The conversion is

    .. math::

        x = r\cos(\phi)\cos(\theta), \qquad
        y = r\cos(\phi)\sin(\theta), \qquad
        z = r\sin(\phi).

    If ``r`` is omitted, a unit sphere is used. The inputs are broadcast by
    NumPy, so scalar and array-valued angles can be combined.

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
    tuple of numpy.ndarray
        The Cartesian coordinates ``(x, y, z)`` with the same shape and
        radial unit as ``r``.
    """
    if r is None:
        r = numpy.ones_like(phi)

    x = r * numpy.cos(phi) * numpy.cos(theta)
    y = r * numpy.cos(phi) * numpy.sin(theta)
    z = r * numpy.sin(phi)
    return x, y, z


def rotate_on_x(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    r"""Rotate spherical coordinates around the X-axis.

    The rotation is applied to the corresponding Cartesian vectors using

    .. math::

        \begin{bmatrix}x'\\y'\\z'\end{bmatrix} =
        \begin{bmatrix}1&0&0\\0&\cos\alpha&-\sin\alpha\\
        0&\sin\alpha&\cos\alpha\end{bmatrix}
        \begin{bmatrix}x\\y\\z\end{bmatrix}.

    The radius is unchanged and is returned together with the rotated angles.

    Parameters
    ----------
    phi : numpy.ndarray
        Azimuthal angles in radians.
    theta : numpy.ndarray
        Polar angles in radians.
    angle : float
        Rotation angle about the X-axis, in radians.

    Returns
    -------
    tuple
        The rotated spherical coordinates ``(r, phi, theta)``.
    """
    # Convert to Cartesian for rotation
    x, y, z = spherical_to_cartesian(phi, theta)
    # Apply rotation around the X-axis
    yp = y * numpy.cos(angle) - z * numpy.sin(angle)
    zp = y * numpy.sin(angle) + z * numpy.cos(angle)
    # Convert back to spherical coordinates
    return cartesian_to_spherical(x, yp, zp)
