
import numpy as np
import fibermodes

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from mpl_toolkits.axes_grid1 import make_axes_locatable

from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Representations import Detector_FarField, LP_FarField
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
                 Filter:            float  = 'None',
                 Name:              str    = 'Intensity Detector'):

        self._name, self._coupling  = Name, 'Intensity'

        self.Source, self.Npts = Source, Npts

        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset

        self.FarField = Detector_FarField(self.Npts, NA, ThetaOffset, PhiOffset)

        self._Filter = Angle(Filter)



    @property
    def ThetaBound(self):
        return self.__ThetaBound

    @ThetaBound.setter
    def ThetaBound(self, val: list):
        self.FarField.ThetaBound = val

    @property
    def Filter(self):
        return self._Filter

    @Filter.setter
    def Filter(self, val):
        self._Filter = Angle(val)

    @property
    def PhiBound(self):
        return self.FarField.__PhiBound

    @PhiBound.setter
    def PhiBound(self, val: list):
        self.FarField.PhiBound = val

    @property
    def PhiOffset(self):
        return self.FarField.PhiOffset

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.FarField.PhiOffset = val

    @property
    def ThetaOffset(self):
        return self.FarField.__ThetaOffset

    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.FarField.ThetaOffset = val

    @property
    def NA(self):
        return self.FarField._NA

    @NA.setter
    def NA(self, val):
        if val >= 1:
            val = 1
        if val <= 0:
            val = 0
        self.FarField.NA = val



    def Coupling(self,
                 Scatterer,
                 Mode         = 'Centered'):

        return Coupling(Scatterer    = Scatterer,
                        Detector     = self,
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

        self.GetFarField()

        if PhiOffset != 0: self.FarField.PhiOffset = PhiOffset

        if ThetaOffset != 0: self.FarField.ThetaOffset = ThetaOffset



    def GetFarField(self):

        Fiber, CoreDiameter = SMF28()

        temp = fibermodes.field.Field(Fiber.source,
                                      fibermodes.Mode(fibermodes.ModeFamily.HE, *self.ModeNumber),
                                      self.Source.Wavelength,
                                      10*CoreDiameter,
                                      self.Npts).Ex()

        temp = np.array(temp, copy=False)

        temp /= (temp.__abs__()).sum()

        if self.Orientation == 'h': temp = temp.T

        self.NearField = _NearField(temp, 10*CoreDiameter, self.Npts)

        temp = np.fft.fft2(self.NearField.Cartesian)

        temp /= self.GenShift()

        temp = np.fft.fftshift(temp)

        temp /= (temp.__abs__()).sum()

        self.FarField = LP_FarField(temp, 10*CoreDiameter, self.Npts, self._NA)


    def GenShift(self):

        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(self.Npts)*(self.Npts-1)/self.Npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T


    @property
    def ThetaBound(self):
        return self.__ThetaBound

    @ThetaBound.setter
    def ThetaBound(self, val: list):
        self.FarField.ThetaBound = val

    @property
    def Filter(self):
        return self._Filter

    @Filter.setter
    def Filter(self, val):
        self._Filter = Angle(val)

    @property
    def PhiBound(self):
        return self.FarField.__PhiBound

    @PhiBound.setter
    def PhiBound(self, val: list):
        self.FarField.PhiBound = val

    @property
    def PhiOffset(self):
        return self.FarField.__PhiOffset

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.FarField.PhiOffset = val

    @property
    def ThetaOffset(self):
        return self.FarField.__ThetaOffset

    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.FarField.ThetaOffset = val

    @property
    def NA(self):
        return self.FarField._NA

    @NA.setter
    def NA(self, val):
        if val >= 1:
            val = 1
        if val <= 0:
            val = 0
        self.FarField.NA = val



    def Coupling(self,
                 Scatterer,
                 Mode         = 'Centered'):

        return Coupling(Scatterer    = Scatterer,
                        Detector     = self,
                        Mode         = Mode)


    def Footprint(self, Scatterer):
        return GetFootprint(Scatterer    = Scatterer,
                            Detector     = self)









def SMF28():
    CoreDiameter = 8.2e-6
    cladDiameter = 125e-6

    Fiber = fiber(core_radius = CoreDiameter,
                  core_index  = 1.4456,
                  clad_radius = cladDiameter,
                  clad_index  = 1.4444)

    return Fiber, CoreDiameter





class _NearField(object):

    def __init__(self, Input, Size, Npts):
        self.Cartesian = Input
        self.Size = Size
        self.Npts = Npts

        self.Meshes = DirectMeshes(Npts   = self.Npts,
                                  XBound = [-self.Size/2, self.Size/2],
                                  YBound = [-self.Size/2, self.Size/2],
                                  XNpts  = self.Npts,
                                  YNpts  = self.Npts)


    def Plot(self):
        fig = plt.figure(figsize=(6,3))
        ax0 = fig.add_subplot(121)
        ax1 = fig.add_subplot(122)

        ax0.pcolormesh(self.Meshes.X.Vector,
                       self.Meshes.Y.Vector,
                       self.Cartesian.real,
                       shading='auto')

        ax0.set_title('Real Part \n Near-Field Cartesian Coordinates')
        ax0.set_xlabel(r'X-Distance x [$\mu$m]')
        ax0.set_ylabel('Y-Distance y  [$\mu$m]')

        ax1.pcolormesh(self.Meshes.X.Vector,
                       self.Meshes.Y.Vector,
                       self.Cartesian.imag,
                       shading='auto')

        ax1.set_title('Imaginary Part \n Near-Field Cartesian Coordinates')
        ax1.set_xlabel(r'X-Distance x  [$\mu$m]')
        ax1.set_ylabel(r'Y-Distance y  [$\mu$m]')
        fig.tight_layout()







class Angle(object):

    def __init__(self, input):

        if input != 'None':
            self.Degree = input
            self.Radian = deg2rad(input)
        else:
            self.Degree = 'None'
            self.Radian = 'None'





# --
