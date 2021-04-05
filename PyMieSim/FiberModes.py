import numpy as np
import fibermodes

from PyMieSim.Physics import FraunhoferDiffraction
from PyMieSim.utils import Normalize


LPList = [(0,1),
          (0,2),
          (0,3),
          (1,1),
          (1,2),
          (1,3),
          (2,1),
          (2,2),
          (3,1),
          (3,2),
          (4,1),
          (5,1)]


def SMF28(mode, Num):
    """Function return an instance of the fiber class specific for a
    SMF28 fiber optic .

    """
    import fibermodes

    CoreDiameter = 8.2e-6
    cladDiameter = 125e-6

    Fiber = fiber()

    SFactor = 100


    Field = fibermodes.field.Field(Fiber.source,
                                  fibermodes.Mode(fibermodes.ModeFamily.HE, mode[0]+1, mode[1]),
                                  940e-9,
                                  Fiber.CoreDiameter*Num/SFactor,
                                  Num).Ex()

    return np.array(Field, copy=False)


def GenLPfiles(LPList, Num=251):
    """Function generate numpy files containing the LP mode field.
    The file directory is: "PyMieSim/LPmodes/LP*.npy"

    Parameters
    ----------
    LPList : :class:`list`
        List of the modes to be computed.
    Num : :class:`int`
        Number of points to evaluate the mode field.

    """
    for mode in LPList:

        filename = f'PyMieSim/LPmodes/LP{mode[0]}{mode[1]}.npy'

        modeField = SMF28(mode = (mode[0], mode[1]) , Num=Num)

        np.save(filename, modeField)

        print(f'Mode LP{mode[0]}{mode[1]} done!')

    print('Files are saved in "PyMieSim/LPmodes/LP*.npy" ')    


def Genfiles(LPList, padWidth = 2000, Num=251):
    """Function generate numpy files containing the FarField of the LP modes.
    The file directory is: "PyMieSim/LPmodes/FLP*.npy"

    Parameters
    ----------
    LPList : :class:`list`
        List of the modes to be computed.
    padWidth : :class:`int`
        The padding for the fourier transform, the higher the larger is the farfield.
    Num : :class:`int`
        Number of points to evaluate the mode field.

    """
    PAD = ((padWidth,padWidth),(padWidth,padWidth))

    for mode in LPList:

        filename = f'PyMieSim/LPmodes/FLP{mode[0]}{mode[1]}.npy'

        FmodeField = SMF28(mode = (mode[0], mode[1]) , Num=Num)

        FmodeField = np.pad(array    = FmodeField,
                            pad_width = PAD,
                            mode      = 'constant')

        FmodeField = FraunhoferDiffraction(FmodeField)[padWidth:-padWidth, padWidth:-padWidth]

        FmodeField = Normalize(FmodeField)

        np.save(filename, FmodeField)

        print(f'Fourier Mode LP{mode[0]}{mode[1]} done!')

    print('Files are saved in "PyMieSim/LPmodes/FLP*.npy" ')


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


if __name__=='__main__':
    GenLPFourierfiles(LPList)
