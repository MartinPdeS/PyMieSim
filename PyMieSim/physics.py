#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


def FraunhoferDiffraction(near_field: numpy.ndarray) -> numpy.ndarray:
    r"""Function compute the Far-Field of a given Near-Field using Fraunhofer
    relation under approximation of small angle.

    Parameters
    ----------
    near_field : :class:`numpy.ndarray`
        Near-Field input [2D].

    Returns
    -------
    :class:`numpy.ndarray`
        Far-Field ouptut [2D].

    """
    temp = numpy.fft.fft2(near_field)

    temp /= GenShift(temp.shape[0])

    return numpy.fft.fftshift(temp)


def GenShift(npts: int) -> numpy.ndarray:
    r"""Function generate a complex shift array that has to be multiplied to FFT in order to obtain
    a phase accurate Fourier transform.

    Parameters
    ----------
    npts : :class:`int`
        Number of point (per dimension) of the FFT array.

    Returns
    -------
    :class:`numpy.ndarray`
        Complex array for FFT.

    """
    if npts % 2 == 1:
        phase_shift = numpy.exp(-complex(0, 1) * numpy.pi * numpy.arange(npts) * (npts - 1) / npts)

        shift_grid, _ = numpy.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T

    else:
        phase_shift = numpy.exp(-complex(0, 1) * numpy.pi * numpy.arange(npts) * (npts) / npts)

        shift_grid, _ = numpy.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T
