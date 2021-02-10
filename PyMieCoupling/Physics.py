import numpy as np
import fibermodes



def FraunhoferDiffraction(nearField):
    """Function compute the Far-Field of a given Near-Field using Fraunhofer
    relation under approximation of small angle.

    Parameters
    ----------
    nearField : np.ndarray
        Near-Field input [2D].

    Returns
    -------
    np.ndarray
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
    npts : int
        Number of point (per dimension) of the FFT array.

    Returns
    -------
    type
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


class _Polarization(object):

    def __init__(self, input):
        if input == None:
            self.Degree = None
            self.Radian = None
        else:
            self.Degree = input
            self.Radian = np.deg2rad(input)




class Source(object):

    def __init__(self,
                 Wavelength:   float,
                 Polarization: float,
                 Power:        float = 1,
                 Radius:       float = 1):

        self.Wavelength = Wavelength

        self.k = 2 * np.pi / Wavelength

        self.Power = Power

        self.Polarization = _Polarization(Polarization)



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


def SMF28():
    """Function return an instance of the fiber class specific for a
    SMF28 fiber optic .

    """
    CoreDiameter = 8.2e-6
    cladDiameter = 125e-6

    Fiber = fiber()

    return Fiber



class fiber(object):
    """Class generating a fiber object from fibermodes package
    (see requirement.txt).

    Parameters
    ----------
    core_radius : float
        Radius of the core of the fiber.
    core_index : float
        Index of the core of the fiber.
    clad_radius : float
        Radius of the clad of the fiber.
    clad_index : float
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
