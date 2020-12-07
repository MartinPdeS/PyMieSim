
import numpy as np
import fibermodes

from PyMieCoupling.functions.converts import Direct2Angle, NA2Angle, Angle2Direct
from PyMieCoupling.classes.Meshes import ScatMeshes
from PyMieCoupling.classes.Fields import Source, LPField, LPFourier
from PyMieCoupling.functions.converts import deg2rad
from PyMieCoupling.functions.Couplings import Coupling, GetFootprint

import fibermodes




class fiber(object):

    def __init__(self,
                 core_radius,
                 core_index,
                 clad_radius,
                 clad_index):

        self.MaxDirect = 2 * clad_radius

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
                 Filter:            float  = 0,
                 Name:              str    = 'Intensity Detector'):

        self._name = Name

        self._coupling = 'Intensity'

        self.NA = NA

        self.Source = Source

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.Npts = Npts

        self.Detection = None

        self.GenMeshes()

        item = np.ones( self.Meshes.Theta.Mesh.Degree.shape )
        item /= (item.shape[0] * item.shape[1])

        self.Fourier = LPFourier(array = item, Meshes = self.Meshes)

        self.Filter, self.FilterRad = Filter, deg2rad(Filter)


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


    def Coupling(self,
                 Scatterer,
                 Polarization = 'NoFiltered',
                 Mode         = 'Centered'):

        return Coupling(Scatterer    = Scatterer,
                        Detector     = self,
                        Polarization = Polarization,
                        Mode         = Mode)


    def Footprint(self, Scatterer):

        return GetFootprint(Scatterer    = Scatterer,
                            Detector     = self)



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
                 NA:            float  = 0.2,
                 Npts:          int   = 101,
                 ThetaOffset:   float = 0,
                 PhiOffset:     float = 0,
                 Filter:        float = 0,
                 Name:          str   = 'Amplitude detector'):

        self._name, self._coupling, self.Fiber = Name, 'Amplitude', Fiber

        self.Filter, self.FilterRad = Filter, deg2rad(Filter)

        Mode = Mode[0]+1, Mode[1]

        self.Mode = fibermodes.Mode(fibermodes.ModeFamily.HE, *Mode)

        self.NA = NA

        FourierScaleFactor = self.Fiber.MaxDirect/Npts/5e-7

        self.__ThetaBound = np.arcsin(self.NA)

        self.Source = Source

        self.Npts = Npts

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.DirectVec = np.linspace(-self.Fiber.MaxDirect/FourierScaleFactor,
                                      self.Fiber.MaxDirect/FourierScaleFactor,
                                      self.Npts)


        self.GenMeshes()

        item = GetLP(Fiber      = self.Fiber.source,
                     Mode       = self.Mode,
                     Wavelength = self.Source.Wavelength,
                     Size       = self.DirectVec[0],
                     Npts       = self.Npts)

        self.Field, self.Fourier = LPField(array = item[0], DirectVec = self.DirectVec), LPFourier(array = item[1], Meshes = self.Meshes)


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


    def Coupling(self,
                 Scatterer,
                 Polarization = 'NoFiltered',
                 Mode         = 'Centered'):

        return Coupling(Scatterer    = Scatterer,
                        Detector     = self,
                        Polarization = Polarization,
                        Mode         = Mode)


    def Footprint(self, Scatterer):

        return GetFootprint(Scatterer    = Scatterer,
                            Detector     = self)




def GetLP(Fiber,
          Mode,
          Wavelength: float,
          Size:       float,
          Npts:       int):

    Field = fibermodes.field.Field(Fiber,
                                   Mode,
                                   Wavelength,
                                   Size,
                                   Npts).Ex()

    Field = np.array(Field, copy=False)

    Field /= (Field.__abs__()).sum()

    Fourier = np.fft.fft2(Field)

    Fourier /= GenShift(Npts)

    Fourier = np.fft.fftshift(Fourier)

    Fourier /= (Fourier.__abs__()).sum()

    return Field, Fourier




def GenShift(Npts):

    phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(Npts)*(Npts-1)/Npts)

    shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

    return shift_grid * shift_grid.T











# --
