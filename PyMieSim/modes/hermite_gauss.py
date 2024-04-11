#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from scipy.special import hermite


def get_mode_field(
        coordinate: numpy.ndarray,
        x_number: int,
        y_number: int,
        wavelength: float = 1.55,
        waist_radius: float = 0.3,
        z: float = 0) -> numpy.ndarray:
    """
    Calculate the Hermite-Gaussian mode field amplitude at given unstructured coordinates,
    normalized so that the L2 norm of the amplitude is 1.

    Parameters:
    coords (array-like): 2xN array of coordinates, where the first row are x-coordinates
                         and the second row are y-coordinates.
    x_number, m (int): Hermite-Gaussian mode indices (n for x, m for y).
    wavelength (float): Wavelength of light in micrometers.
    waist_radius (float): Waist radius of the beam at the focus in micrometers.
    z (float): Longitudinal position from the beam waist in micrometers.

    Returns:
    np.ndarray: Array of complex field amplitudes at the given coordinates, normalized to L2 norm of 1.
    """
    k = 2 * numpy.pi / wavelength  # Wave number in vacuum
    w0 = waist_radius  # Beam waist

    # Extract x and y coordinates
    x = coordinate[0, :]
    y = coordinate[1, :]

    # Calculate the beam width at distance z
    w = w0 * numpy.sqrt(1 + (z * wavelength / (numpy.pi * w0**2))**2)

    # Calculate the radius of curvature of the beam's wavefront at z
    R = float('inf') if z == 0 else z * (1 + (numpy.pi * w0**2 / (z * wavelength))**2)

    # Calculate the Gouy phase shift at z
    gouy_phase = numpy.arctan(z * numpy.pi / (wavelength * w0**2))

    # Hermite polynomial factors
    Hn = hermite(x_number)(numpy.sqrt(2) * x / w)
    Hm = hermite(y_number)(numpy.sqrt(2) * y / w)

    # Amplitude calculation
    amplitude = Hn * Hm * numpy.exp(-((x**2 + y**2) / w**2))

    # Phase calculation including Gouy phase and spherical phase factor
    phase = -k * ((x**2 + y**2) / (2 * R)) + (x_number + y_number) * gouy_phase
    field = amplitude * numpy.exp(1j * phase)

    # Normalization to L2 norm of 1
    norm = numpy.sqrt(numpy.sum(numpy.abs(field)**2))
    field /= norm

    return field


def interpolate_from_fibonacci_mesh(fibonacci_mesh, **kwargs) -> numpy.ndarray:
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
        np.ndarray: Array of complex field amplitudes interpolated at the coordinates defined by the Fibonacci mesh.
    """
    coordinate = numpy.row_stack((
        fibonacci_mesh.base_x,
        fibonacci_mesh.base_y,
    ))

    norm = numpy.sqrt(numpy.square(coordinate).sum(axis=0)).max()

    mode_field = get_mode_field(coordinate[:2, :] / norm, **kwargs)

    return mode_field


def interpolate_from_structured_mesh(sampling: int = 50, **kwargs) -> numpy.ndarray:
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

    mode_field = get_mode_field(coordinate[:2, :] / norm, **kwargs)

    return mode_field.reshape([sampling, sampling])

# -
