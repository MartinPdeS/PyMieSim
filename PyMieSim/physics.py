#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from scipy.constants import epsilon_0, c


def fraunhofer_diffraction(near_field: numpy.ndarray) -> numpy.ndarray:
    """
    Calculate the Far-Field diffraction pattern from a given Near-Field pattern using the Fraunhofer approximation.

    Parameters:
        - near_field: numpy.ndarray, the near field pattern as a 2D numpy array.

    Returns:
        - numpy.ndarray, the far field diffraction pattern as a 2D numpy array.
    """
    fft_result = numpy.fft.fft2(near_field)
    fft_result /= generate_fft_shift(fft_result.shape[0])
    far_field = numpy.fft.fftshift(fft_result)
    return far_field


def generate_fft_shift(n_points: int) -> numpy.ndarray:
    """
    Generate a complex shift array for phase correction in FFT, ensuring accurate Fourier transform phase.

    Parameters:
        - n_points: int, number of points in one dimension of the square array.

    Returns:
        - numpy.ndarray, a complex shift array for FFT phase correction.
    """
    phase_arg = numpy.pi * numpy.arange(n_points) * (n_points - (n_points % 2)) / n_points
    phase_shift = numpy.exp(-1j * phase_arg)
    shift_grid = numpy.outer(phase_shift, phase_shift)
    return shift_grid


def power_to_amplitude(optical_power: float, NA: float, wavelength: float) -> float:
    """
    Calculate the optical field amplitude from the given optical power, numerical aperture (NA), and wavelength.

    Parameters:
        - optical_power: float, the optical power in watts.
        - NA: float, the numerical aperture of the optical system.
        - wavelength: float, the wavelength of the light in meters.

    Returns:
        - float, the amplitude of the optical field.
    """
    omega = 0.61 * wavelength / NA
    area = numpy.pi * (omega / 2)**2
    intensity = optical_power / area
    amplitude = numpy.sqrt(2 * intensity / (c * epsilon_0))
    return amplitude
