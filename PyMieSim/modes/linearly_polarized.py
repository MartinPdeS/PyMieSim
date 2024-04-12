#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy
from scipy.special import jn, jn_zeros


def get_mode_field(
        coords,
        azimuthal_number: int,
        radial_number: int,
        core_radius: float = 1,
        wavelength=1.55,
        n_core=1.45,
        n_cladding=1.44):
    """
    Calculate the LP mode field amplitude at given structured coordinates,
    normalized so that the integral of the absolute square of the amplitude over
    the area is 1.

    Parameters:
        coordinates (array-like): 2xN array of coordinates in cylindrical system (r, phi),
            where the first row are radial coordinates, and the second row
            are angular coordinates.
        azimuthal_number (int): Azimuthal index of the LP mode.
        radial_number (int): Radial index of the LP mode.
        core_radius (float): Core radius of the fiber in micrometers.
        wavelength (float): Wavelength of light in micrometers.
        n_core (float): Refractive index of the fiber core.
        n_cladding (float): Refractive index of the fiber cladding.

    Returns:
    numpy.ndarray: Array of complex field amplitudes at the given coordinates, normalized.
    """
    k = 2 * numpy.pi / wavelength  # noqa: F841

    # Bessel function roots
    u = jn_zeros(azimuthal_number, radial_number)[-1]  # Take the m-th zero of the l-th order Bessel function

    # Extract r and phi coordinates
    r = numpy.sqrt(numpy.square(coords).sum(axis=0))
    phi = numpy.arctan2(coords[0], coords[1])

    # Radial part using Bessel function
    radial_part = jn(azimuthal_number, u * r / core_radius)

    # Azimuthal part
    azimuthal_part = numpy.cos(azimuthal_number * phi)

    # Field
    field = radial_part * azimuthal_part

    # Normalization to power of 1
    norm = numpy.sqrt(numpy.sum(numpy.abs(field)**2))
    field /= norm

    return field


def interpolate_from_structured_mesh(number_0: int, number_1: int, sampling: int = 50) -> numpy.ndarray:
    """
    Generate a structured mesh grid.

    Parameters:
        sampling (int): Number of points in each dimension of the grid.

    Returns:
        numpy.ndarray: 2xN array of mesh coordinates [x, y] from -100 to 100.
    """
    x_mesh, y_mesh = numpy.mgrid[-100:100:complex(sampling), -100:100:complex(sampling)]

    coordinate = numpy.row_stack((
        x_mesh.ravel(),
        y_mesh.ravel(),
    ))

    norm = numpy.sqrt(numpy.square(coordinate).sum(axis=0)).max()

    coordinate /= norm

    mode_field = get_mode_field(
        coordinate[:2, :],
        azimuthal_number=number_0,
        radial_number=number_1
    )

    return mode_field.reshape([sampling, sampling])


def interpolate_from_fibonacci_mesh(fibonacci_mesh, number_0: int, number_1: int) -> numpy.ndarray:
    """
    Calculate the Hermite-Gaussian mode field for given mode indices on a Fibonacci mesh.

    Parameters:
        fibonacci_mesh (object): An object with attributes 'x' and 'y' containing the mesh coordinates.
        n (int): Hermite-Gaussian mode index along the x-direction.
        m (int): Hermite-Gaussian mode index along the y-direction.
        wavelength (float): Wavelength of the light in micrometers. Default is 1.55 micrometers.
        waist_radius (float): Waist radius of the beam at the focus in micrometers. Default is 1.0 micrometers.
        z (float): Axial position from the beam waist in micrometers where the field is calculated. Default is 0.

    Returns:
        numpy.ndarray: Array of complex field amplitudes interpolated at the coordinates defined by the Fibonacci mesh.
    """
    coordinate = numpy.row_stack((
        fibonacci_mesh.base_x,
        fibonacci_mesh.base_y,
    ))

    norm = numpy.sqrt(numpy.square(coordinate).sum(axis=0)).max()

    coordinate /= norm

    mode_field = get_mode_field(
        coordinate[:2, :],
        azimuthal_number=number_0,
        radial_number=number_1
    )

    return mode_field
