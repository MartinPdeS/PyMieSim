
import numpy as np
import cupy as cp
import fibermodes

from PyMieCoupling.functions.converts import rad2deg, deg2rad, Angle2Direct, Direct2Angle, NA2Angle
from PyMieCoupling.classes.Meshes import Meshes
from PyMieCoupling.classes.Misc import Source, LPField, LPFourier
from PyMieCoupling.functions.Misc import GetLP


class Photodiode(object):
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
                 GPU:               bool   = False,
                 Name:              str    = 'Detector'):

        self._name = Name

        self.GPU = GPU

        self._coupling = 'Intensity'

        self.NA = NA

        self.Source = Source

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.__ThetaBound = np.arcsin(self.NA)

        self.Npts = Npts

        self.GenMeshes()

        item = self.GenField()

        self.Fourier = LPFourier(item, self.Meshes)


    def GenMeshes(self):
        self.__ThetaBound, self.__PhiBound  = NA2Angle(self.NA, self.GPU)

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = (self.__ThetaBound) + self.__ThetaOffset,
                             PhiBound   = (self.__PhiBound) + self.__PhiOffset,
                             GPU        = self.GPU)


    def GenField(self):
        if self.GPU:
            return cp.ones( cp.shape( self.Meshes.Theta.Mesh.Degree ) )
        else:
            return np.ones( np.shape( self.Meshes.Theta.Mesh.Degree ) )


    @property
    def ThetaBound(self):
        return self.__ThetaBound


    @ThetaBound.setter
    def ThetaBound(self, val: list):

        self.__ThetaBound = val

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             GPU        = self.GPU)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiBound(self):
        return self.__PhiBound


    @PhiBound.setter
    def PhiBound(self, val: list):

        self.__PhiBound = val

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             GPU        = self.GPU)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiOffset(self):
        return self.__PhiOffset


    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val

        self.PhiBound += val



    @property
    def ThetaOffset(self):
        return self.__ThetaOffset


    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.__ThetaOffset = val

        self.ThetaBound += val



class LPmode(object):
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
                 Fiber:         fibermodes.fiber,
                 Mode:          tuple,
                 Source:        Source,
                 Npts:          int   = 101,
                 Magnification: float = 1.,
                 ThetaOffset:   float = 0,
                 PhiOffset:     float = 0,
                 GPU:           bool  = False,
                 Name:          str   = 'Field detector'):

        self._name, self._coupling, self.Fiber = Name, 'Amplitude', Fiber

        self.GPU = GPU

        Mode = Mode[0]+1, Mode[1]

        self.Mode = fibermodes.Mode(fibermodes.ModeFamily.HE, *Mode)

        self.Source = Source

        self.Npts, self.Magnification  = Npts, Magnification

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.DirectVec = np.linspace(-self.Fiber.MaxDirect, self.Fiber.MaxDirect, self.Npts)

        self.GenMeshes()

        item = self.GenField()

        self.Field, self.Fourier = LPField(item[0], self.DirectVec), LPFourier(item[1], self.Meshes)

        if Magnification != 1:
            self.Magnificate(Magnification)


    def GenField(self):

        return GetLP(Fiber      = self.Fiber.source,
                     Mode       = self.Mode,
                     Wavelength = self.Source.Wavelength,
                     Size       = self.DirectVec[0],
                     Npts       = self.Npts,
                     GPU        = self.GPU)


    def GenMeshes(self):

        self.AngleVec = Direct2Angle(self.DirectVec, self.Source.k)

        self._DirectBound = [self.DirectVec[0], self.DirectVec[-1]]

        self.__ThetaBound = np.array( [self.AngleVec[0], self.AngleVec[-1]] )

        self.__PhiBound = np.array( [self.AngleVec[0], self.AngleVec[-1]] )

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound + self.__ThetaOffset,
                             PhiBound   = self.__PhiBound + self.__PhiOffset,
                             GPU        = self.GPU)



    def Magnificate(self, Magnification):

        self.DirectVec /= Magnification

        self.GenMeshes()

        self.Field.DirectVec, self.Fourier.Meshes = self.DirectVec, self.Meshes


    @property
    def ThetaBound(self):
        return self.__ThetaBound


    @ThetaBound.setter
    def ThetaBound(self, val: list):

        self.__ThetaBound = val

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             GPU        = self.GPU)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiBound(self):
        return self.__PhiBound


    @PhiBound.setter
    def PhiBound(self, val: list):

        self.__PhiBound = val

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             GPU        = self.GPU)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiOffset(self):
        return self.__PhiOffset

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val

        self.PhiBound += val


    @property
    def ThetaOffset(self):
        return self.__ThetaOffset

    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.__ThetaOffset = val

        self.ThetaBound += val











# --
