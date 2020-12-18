
import numpy as np


import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from mpl_toolkits.axes_grid1 import make_axes_locatable

from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes
from PyMieCoupling.utils import Source, SMF28, Angle
from PyMieCoupling.classes.Fields import Detector_FarField, LPFarField, LPNearField
from PyMieCoupling.functions.converts import deg2rad
from PyMieCoupling.functions.Couplings import Coupling, GetFootprint
from PyMieCoupling.classes.BaseClasses import BaseDetector
import fibermodes






class Photodiode(BaseDetector):
    """Short summary.

    Parameters
    ----------
    size : float
        Size of the detector, [diameter for circle shaped/ side for square].
    shape : str
        Shape of the detector.
    wavelength : float
        Wavelength of the incoming source.
    npts : int
        Number of points defining the rastered meshes.
    ThetaOffset : float
        Offset of theta angle between the detector and source.
    PhiOffset : float
        Offset of phi angle between the detector and source.
    Magnification : float
        Magnification induced by the lense.
    Name : str
        Name of detector [optional for plots].

    Attributes
    ----------
    _name : type
        Description of attribute `_name`.
    _coupling : type
        Description of attribute `_coupling`.
    k : type
        Description of attribute `k`.
    DirectVec : type
        Description of attribute `DirectVec`.
    GenShift : type
        Description of attribute `GenShift`.
    GenMeshes : type
        Description of attribute `GenMeshes`.
    Field : type
        Description of attribute `Field`.
    Fourier : type
        Description of attribute `Fourier`.
    GenField : type
        Description of attribute `GenField`.
    magnificate : type
        Description of attribute `magnificate`.
    size
    wavelength
    ThetaOffset
    PhiOffset
    npts

    """

    def __init__(self,
                 Source:            Source = None,
                 NA:                float  = 0.2,
                 Npts:              int    = 101,
                 ThetaOffset:       float  = 0,
                 PhiOffset:         float  = 0,
                 Filter:            float  = 'None',
                 Name:              str    = 'Intensity Detector'):

        self._name, self._coupling  = Name, 'Intensity'

        self.Source, self.Npts = Source, Npts

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.FarField = Detector_FarField(self.Npts, NA, ThetaOffset, PhiOffset)

        self._Filter = Angle(Filter)





class LPmode(BaseDetector):
    """Short summary.

    Parameters
    ----------
    Fiber : float
        The Fiber class used to generate the modes.
    LPmode : str
        The indices of LP mode.
    wavelength : float
        Wavelength of the incoming source.
    npts : int
        Number of points defining the rastered meshes.
    ThetaOffset : float
        Offset of theta angle between the detector and source.
    PhiOffset : float
        Offset of phi angle between the detector and source.
    Magnification : float
        Magnification induced by the lense.
    Name : str
        Name of detector [optional for plots].

    Attributes
    ----------
    _name : type
        Description of attribute `_name`.
    _coupling : type
        Description of attribute `_coupling`.
    k : type
        Description of attribute `k`.
    DirectVec : type
        Description of attribute `DirectVec`.
    GenMeshes : type
        Description of attribute `GenMeshes`.
    Field : type
        Description of attribute `Field`.
    Fourier : type
        Description of attribute `Fourier`.
    GenField : type
        Description of attribute `GenField`.
    magnificate : type
        Description of attribute `magnificate`.
    size
    wavelength
    ThetaOffset
    PhiOffset
    npts

    """
    def __init__(self,
                 Mode:          tuple,
                 Source:        Source,
                 Orientation:   str   = 'h',
                 NA:            float = 0.2,
                 Npts:          int   = 101,
                 ThetaOffset:   float = 0,
                 PhiOffset:     float = 0,
                 Filter:        float =  'None',
                 Name:          str   = 'Amplitude detector'):

        assert Orientation in ['v','h'], "Orientation should either be v [vertical] or h [horizontal]'"

        self._name, self._coupling = Name, 'Amplitude'

        self._NA = NA

        self._Filter = Angle(Filter)

        self.ModeNumber = Mode[0]+1, Mode[1]

        self.Source, self.Npts = Source, Npts

        self.Orientation = Orientation

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        Fiber, CoreDiameter = SMF28()

        if PhiOffset == 0:

            self.GetFarFieldNoOffset()



    def GetFarFieldNoOffset(self):

        Fiber, CoreDiameter = SMF28()

        temp = fibermodes.field.Field(Fiber.source,
                                      fibermodes.Mode(fibermodes.ModeFamily.HE, *self.ModeNumber),
                                      self.Source.Wavelength,
                                      10*CoreDiameter,
                                      self.Npts).Ex()

        temp = np.array(temp, copy=False)

        temp /= (temp.__abs__()).sum()

        if self.Orientation == 'h': temp = temp.T

        self.NearField = LPNearField(temp, 10*CoreDiameter, self.Npts)

        temp = np.fft.fft2(self.NearField.Cartesian)

        temp /= self.GenShift()

        temp = np.fft.fftshift(temp)

        temp /= (temp.__abs__()).sum()

        self.FarField = LPFarField(temp, 10*CoreDiameter, self.Npts, self._NA)


    def GenShift(self):

        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(self.Npts)*(self.Npts-1)/self.Npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T















# --
