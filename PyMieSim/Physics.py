#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import fibermodes

class _Polarization(object):

    def __init__(self, input):
        if input == None:
            self.Degree = None
            self.Radian = None
        else:
            self.Degree = input
            self.Radian = np.deg2rad(input)

def FraunhoferDiffraction(nearField):
    """Function compute the Far-Field of a given Near-Field using Fraunhofer
    relation under approximation of small angle.

    Parameters
    ----------
    nearField : :class:`np.ndarray`
        Near-Field input [2D].

    Returns
    -------
    :class:`np.ndarray`
        Far-Field ouptut [2D].

    """

    temp = np.fft.fft2(nearField)

    temp /= GenShift(temp.shape[0])

    return np.fft.fftshift(temp)



def GenShift(npts):
    """Function generate a complex shift array that has to be multiplied to FFT in order to obtain
    a phase accurate Fourier transform.

    Parameters
    ----------
    npts : :class:`int`
        Number of point (per dimension) of the FFT array.

    Returns
    -------
    :class:`np.ndarray`
        Complex array for FFT.

    """

    if npts % 2 == 1 :
        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(npts)*(npts-1)/npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T

    else:
        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(npts)*(npts)/npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T



def SMF28(mode, Num, SFactor):
    """Function return an instance of the fiber class specific for a
    SMF28 fiber optic .

    """
    CoreDiameter = 8.2e-6
    cladDiameter = 125e-6

    Fiber = fiber()

    Field = fibermodes.field.Field(Fiber.source,
                                  fibermodes.Mode(fibermodes.ModeFamily.HE, mode[0]+1, mode[1]),
                                  940e-9,
                                  Fiber.CoreDiameter*Num/SFactor,
                                  Num).Ex()

    return np.array(Field, copy=False)


class Angle(object):

    def __init__(self, input, unit='Degree'):
        if np.asarray(input).any() == None:
            self.Degree = None
            self.Radian = None

        if unit == 'Degree':
            self.Degree = input
            self.Radian = np.deg2rad(input)
        if unit == 'Radian':
            self.Degree = np.rad2deg(input)
            self.Radian = input


class fiber(object):
    """Class generating a fiber object from fibermodes package
    (see requirement.txt).

    Parameters
    ----------
    core_radius : :class:`float`
        Radius of the core of the fiber.
    core_index : :class:`float`
        Index of the core of the fiber.
    clad_radius : :class:`float`
        Radius of the clad of the fiber.
    clad_index : :class:`float`
        Index of the clad of the fiber.

    """

    def __init__(self,
                 core_radius: float  = 8.2e-6,
                 core_index:  float  = 1.4456,
                 clad_radius: float  = 125e-6,
                 clad_index:  float  = 1.4444):

        self.MaxDirect = 2 * clad_radius

        self.CoreDiameter = core_radius

        factory = fibermodes.FiberFactory()

        factory.addLayer(name     = 'core',
                         radius   = core_radius,
                         material = 'Fixed',
                         geometry = "StepIndex",
                         index    = 1.4489)

        factory.addLayer(name     = 'cladding',
                         material = 'Fixed',
                         index    = 1)

        self.source = factory[0]
