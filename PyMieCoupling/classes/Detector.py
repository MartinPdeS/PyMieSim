
import numpy as np
import fibermodes

from PyMieCoupling.functions.converts import Direct2Angle, NA2Angle
from PyMieCoupling.classes.Meshes import ScatMeshes
from PyMieCoupling.functions.Misc import GetLP
from PyMieCoupling.classes.Misc import Source, LPField, LPFourier


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
                 filter:            float  = None,
                 Name:              str    = 'Detector'):

        self._name = Name

        self._coupling = 'Intensity'

        self.NA = NA

        self.Source = Source

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.__ThetaBound = np.arcsin(self.NA)

        self.Npts = Npts

        self.Detection = None

        self.GenMeshes()

        item = np.ones( self.Meshes.Theta.Mesh.Degree.shape ) / (self.Meshes.Theta.Mesh.Degree.shape[0]*self.Meshes.Theta.Mesh.Degree.shape[1])

        self.Fourier = LPFourier(array = item, Meshes = self.Meshes)

        self.filter = filter


    def GenMeshes(self):
        self.__ThetaBound, self.__PhiBound  = NA2Angle(self.NA)

        self.Meshes = ScatMeshes(Npts       = self.Npts,
                                 ThetaBound = (self.__ThetaBound) + self.__ThetaOffset,
                                 PhiBound   = (self.__PhiBound) + self.__PhiOffset)



    @property
    def ThetaBound(self):
        return self.__ThetaBound


    @ThetaBound.setter
    def ThetaBound(self, val: list):

        self.__ThetaBound = val

        self.Meshes = ScatMeshes(Npts       = self.Npts,
                                 ThetaBound = self.__ThetaBound,
                                 PhiBound   = self.__PhiBound)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiBound(self):
        return self.__PhiBound


    @PhiBound.setter
    def PhiBound(self, val: list):

        self.__PhiBound = val

        self.Meshes = ScatMeshes(Npts       = self.Npts,
                                 ThetaBound = self.__ThetaBound,
                                 PhiBound   = self.__PhiBound)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiOffset(self):
        return self.__PhiOffset


    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val

        self.PhiBound = np.array([-self.Meshes.Phi.Range.Degree/2, self.Meshes.Phi.Range.Degree/2], copy=False) + val



    @property
    def ThetaOffset(self):
        return self.__ThetaOffset


    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.__ThetaOffset = val

        self.ThetaBound = np.array([-self.Meshes.Theta.Range.Degree/2, self.Meshes.Theta.Range.Degree/2], copy=False) + val


    def Coupling(self, Source):


        dOmega = self.Meshes.Phi.Delta.Radian *\
                 self.Meshes.Theta.Delta.Radian

        Perp = self.Fourier *\
               (Source.Field.Perpendicular).__abs__() *\
               (np.sin(self.Meshes.Phi.Mesh.Radian + np.pi/2).T).__abs__()


        Para = self.Fourier *\
               (Source.Field.Parallel).__abs__() *\
               (np.sin(self.Meshes.Phi.Mesh.Radian + np.pi/2).T).__abs__()

        if self.filter:
            PerpFiltre = (np.cos(Perp) * dOmega).sum()**2
            ParaFiltre = (np.sin(Para) * dOmega).sum()**2
            Detection = PerpFiltre + ParaFiltre

        else:
            PerpFiltre = (Perp * dOmega).sum()**2
            ParaFiltre = (Para * dOmega).sum()**2
            Detection = PerpFiltre + ParaFiltre

        Perp = (Perp * dOmega).sum()**2

        Para = (Para * dOmega).sum()**2

        return {'Parallel': Para, 'Perpendicular': Perp, 'Filtered': Detection}




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
                 filter:        float = None,
                 Name:          str   = 'Field detector'):

        self._name, self._coupling, self.Fiber = Name, 'Amplitude', Fiber

        self.Filter = filter

        Mode = Mode[0]+1, Mode[1]

        self.Mode = fibermodes.Mode(fibermodes.ModeFamily.HE, *Mode)

        self.Source = Source

        self.Npts, self.Magnification  = Npts, Magnification

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.DirectVec = np.linspace(-self.Fiber.MaxDirect, self.Fiber.MaxDirect, self.Npts)

        self.GenMeshes()

        item = GetLP(Fiber      = self.Fiber.source,
                     Mode       = self.Mode,
                     Wavelength = self.Source.Wavelength,
                     Size       = self.DirectVec[0],
                     Npts       = self.Npts)

        self.Field, self.Fourier = LPField(array = item[0], DirectVec = self.DirectVec), LPFourier(array = item[1], Meshes = self.Meshes)

        if Magnification != 1:
            self.Magnificate(Magnification)


    def GenMeshes(self):

        self.AngleVec = np.array( Direct2Angle(self.DirectVec, self.Source.k), copy=False )

        self._DirectBound = [self.DirectVec[0], self.DirectVec[-1]]

        self.__ThetaBound = self.AngleVec + self.__PhiOffset#

        self.__ThetaBound = np.array( [ self.AngleVec[0], self.AngleVec[-1] ], copy=False ) + self.__ThetaOffset

        self.__PhiBound = np.array( [ self.AngleVec[0], self.AngleVec[-1] ], copy=False ) + self.__PhiOffset

        self.Meshes = ScatMeshes(Npts       = self.Npts,
                                 ThetaBound = self.__ThetaBound,
                                 PhiBound   = self.__PhiBound)



    def Magnificate(self, Magnification):

        self.DirectVec /= Magnification

        self.GenMeshes()

        self.Field.DirectVec, self.Fourier.Meshes = self.DirectVec, self.Meshes


    @property
    def ThetaBound(self):
        return self.__ThetaBound


    @ThetaBound.setter
    def ThetaBound(self, val: list):

        self.__ThetaBound = np.array( val, copy=False )

        self.Meshes = ScatMeshes(Npts       = self.Npts,
                                 ThetaBound = self.__ThetaBound,
                                 PhiBound   = self.__PhiBound)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiBound(self):
        return self.__PhiBound


    @PhiBound.setter
    def PhiBound(self, val: list):

        self.__PhiBound = np.array( val, copy=False )

        self.Meshes = ScatMeshes(Npts       = self.Npts,
                                 ThetaBound = self.__ThetaBound,
                                 PhiBound   = self.__PhiBound)

        self.Fourier.Meshes = self.Meshes


    @property
    def PhiOffset(self):
        return self.__PhiOffset


    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val
        self.PhiBound = np.array([-self.Meshes.Theta.Range.Degree/2, self.Meshes.Theta.Range.Degree], copy=False) + val


    @property
    def ThetaOffset(self):
        return self.__ThetaOffset


    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.__ThetaOffset = val
        self.ThetaBound = np.array([-self.Meshes.Theta.Range.Degree/2, self.Meshes.Theta.Range.Degree/2], copy=False) + val


    def Coupling(self, Source):

        dOmega = self.Meshes.Phi.Delta.Radian *\
                 self.Meshes.Theta.Delta.Radian

        Perp = self.Field *\
               Source.Field.Perpendicular *\
               np.sin(self.Meshes.Phi.Mesh.Radian.T).__abs__() *\
               dOmega

        CPerp = Perp.sum().__abs__()**2

        Para = self.Field *\
               Source.Field.Parallel *\
               np.sin(self.Meshes.Phi.Mesh.Radian.T).__abs__() *\
               dOmega

        CPara = Para.sum().__abs__()**2


        if self.Filter:
            Filtered = (Perp * np.sin(self.Filter/180*np.pi) + Para * np.cos(self.Filter/180*np.pi)).sum().__abs__()**2
        else:
            Filtered = (Perp + Para).sum().__abs__()**2


        return {'Parallel': CPara, 'Perpendicular': CPerp, 'Filtered': Filtered}



























# --
