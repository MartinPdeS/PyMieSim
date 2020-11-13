
import numpy as np
import fibermodes

from PyMieCoupling.functions.converts import rad2deg, deg2rad, Angle2Direct, Direct2Angle, NA2Angle
from PyMieCoupling.classes.Meshes import Meshes
from PyMieCoupling.classes.Misc import Source, LPField, LPFourier, Operation as Op
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
                 cuda:              bool   = False,
                 Name:              str    = 'Detector'):

        self._name = Name

        self.cuda = cuda

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
        self.__ThetaBound, self.__PhiBound  = NA2Angle(self.NA, self.cuda)

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = (self.__ThetaBound) + self.__ThetaOffset,
                             PhiBound   = (self.__PhiBound) + self.__PhiOffset,
                             cuda       = self.cuda)


    def GenField(self):
            return Op.ones(self.cuda)( self.Meshes.Theta.Mesh.Degree.shape )


    @property
    def ThetaBound(self):
        return self.__ThetaBound


    @ThetaBound.setter
    def ThetaBound(self, val: list):

        self.__ThetaBound = val

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             cuda       = self.cuda)

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
                             cuda       = self.cuda)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiOffset(self):
        return self.__PhiOffset


    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val

        self.PhiBound = np.array([-self.Meshes.Phi.Range.Degree/2, self.Meshes.Phi.Range.Degree/2]) + val



    @property
    def ThetaOffset(self):
        return self.__ThetaOffset


    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.__ThetaOffset = val

        self.ThetaBound = np.array([-self.Meshes.Theta.Range.Degree/2, self.Meshes.Theta.Range.Degree/2]) + val


    def Coupling(self, Source):

        dOmega = self.Meshes.Phi.Delta.Radian *\
                 self.Meshes.Theta.Delta.Radian

        Perp = self.Fourier.Array *\
               (Source.Field.Perpendicular).__abs__() *\
               (Op.sin(self.Meshes.Phi.Mesh.Radian + Op.pi/2).T).__abs__()

        Perp = (Perp * dOmega).sum()**2

        Para = self.Fourier.Array *\
               (Source.Field.Parallel).__abs__() *\
               (Op.sin(self.Meshes.Phi.Mesh.Radian + Op.pi/2).T).__abs__()

        Para = (Para * dOmega).sum()**2

        return {'Parallel': Para, 'Perpendicular': Perp}




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
                 cuda:          bool  = False,
                 Name:          str   = 'Field detector'):

        self._name, self._coupling, self.Fiber, self.cuda = Name, 'Amplitude', Fiber, cuda

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
                     cuda       = self.cuda)


    def GenMeshes(self):

        self.AngleVec = np.array( Direct2Angle(self.DirectVec, self.Source.k, cuda = False) )

        self._DirectBound = [self.DirectVec[0], self.DirectVec[-1]]

        self.__ThetaBound = self.AngleVec + self.__PhiOffset#

        self.__ThetaBound = np.array( [ self.AngleVec[0], self.AngleVec[-1] ] ) + self.__ThetaOffset

        self.__PhiBound = np.array( [ self.AngleVec[0], self.AngleVec[-1] ] ) + self.__PhiOffset

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             cuda       = self.cuda)



    def Magnificate(self, Magnification):

        self.DirectVec /= Magnification

        self.GenMeshes()

        self.Field.DirectVec, self.Fourier.Meshes = self.DirectVec, self.Meshes


    @property
    def ThetaBound(self):
        return self.__ThetaBound


    @ThetaBound.setter
    def ThetaBound(self, val: list):

        self.__ThetaBound = np.array( val )

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             cuda       = self.cuda)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiBound(self):
        return self.__PhiBound


    @PhiBound.setter
    def PhiBound(self, val: list):

        self.__PhiBound = np.array( val )

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self.__ThetaBound,
                             PhiBound   = self.__PhiBound,
                             cuda       = self.cuda)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiOffset(self):
        return self.__PhiOffset


    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val
        self.PhiBound = np.array([-self.Meshes.Theta.Range.Degree/2, self.Meshes.Theta.Range.Degree]) + val


    @property
    def ThetaOffset(self):
        return self.__ThetaOffset


    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.__ThetaOffset = val
        self.ThetaBound = np.array([-self.Meshes.Theta.Range.Degree/2, self.Meshes.Theta.Range.Degree/2]) + val



    def Coupling(self, Source):

        dOmega = self.Meshes.Phi.Delta.Radian *\
                 self.Meshes.Theta.Delta.Radian

        Perp = self.Field.Array *\
               Source.Field.Perpendicular *\
               (Op.sin(self.Meshes.Phi.Mesh.Radian).T).__abs__()

        Perp = ( Perp.sum() ).__abs__()**2

        Para = self.Field.Array *\
               Source.Field.Parallel *\
               (Op.sin(self.Meshes.Phi.Mesh.Radian).T).__abs__()

        Para = (Para * dOmega).sum().__abs__()**2

        return {'Parallel': Para, 'Perpendicular': Perp}






# --
