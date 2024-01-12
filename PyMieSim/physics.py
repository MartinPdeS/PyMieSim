#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


def FraunhoferDiffraction(near_field: numpy.ndarray) -> numpy.ndarray:
    """
    Function compute the Far-Field of a given Near-Field using Fraunhofer
    relation under approximation of small angle.

    :param      near_field:  The near field
    :type       near_field:  numpy.ndarray

    :returns:   the far field
    :rtype:     numpy.ndarray
    """
    temp = numpy.fft.fft2(near_field)

    temp /= generate_fft_shift(temp.shape[0])

    return numpy.fft.fftshift(temp)


def generate_fft_shift(n_points: int) -> numpy.ndarray:
    """
    Function generate a complex shift array that has to be multiplied to FFT in order to obtain
    a phase accurate Fourier transform.

    :param      npts:  The npts
    :type       npts:  int

    :returns:   Complex array for FFT.
    :rtype:     numpy.ndarray
    """
    if n_points % 2 == 1:
        phase_shift = numpy.exp(-complex(0, 1) * numpy.pi * numpy.arange(n_points) * (n_points - 1) / n_points)

        shift_grid, _ = numpy.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T

    else:
        phase_shift = numpy.exp(-complex(0, 1) * numpy.pi * numpy.arange(n_points) * (n_points) / n_points)

        shift_grid, _ = numpy.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T


def power_to_amplitude(optical_power: float, NA: float, wavelength: float) -> float:
    """
    Returns the amplitude of field for a given power, numerical aperture and wavelength

    :param      power:       The power
    :type       power:       float
    :param      NA:          { parameter_description }
    :type       NA:          float
    :param      wavelength:  The wavelength
    :type       wavelength:  float

    :returns:   { description_of_the_return_value }
    :rtype:     float
    """
    from scipy.constants import epsilon_0, c
    omega = 0.61 * wavelength / NA
    area = numpy.pi * (omega / 2)**2
    intensity = optical_power / area
    amplitude = numpy.sqrt(
        2 / (c * epsilon_0) * intensity
    )

    return amplitude

# -
