
import numpy as np


import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from mpl_toolkits.axes_grid1 import make_axes_locatable


from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes
from PyMieCoupling.classes.Representations import ScalarFarField
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Source, SMF28, Angle, _Polarization
from PyMieCoupling.classes.BaseClasses import BaseFarField
from PyMieCoupling.functions.converts import deg2rad
from PyMieCoupling.classes.BaseClasses import BaseDetector
import fibermodes
from scipy.interpolate import griddata
import polarTransform





class Photodiode(BaseDetector):

    def __init__(self,
                 Source:            Source = None,
                 NA:                float  = 0.2,
                 Samples:           int    = 401,
                 GammaOffset:       float  = 0,
                 PhiOffset:         float  = 0,
                 Filter:            float  = 'None',
                 Name:              str    = 'Intensity Detector'):

        self._name, self._coupling  = Name, 'Intensity'

        self._NA = NA

        self.Source, self.Samples = Source, Samples

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.GetSpherical()


    def GetSpherical(self):


        self.Meshes = AngleMeshes(MaxAngle    = NA2Angle(self._NA)*np.pi/180,
                                  Samples     = self.Samples,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)

        self.Samples = self.Meshes.Phi.Radian.shape

        Scalar = np.ones(self.Samples) / (self.Samples)

        self.Scalar = ScalarFarField(Scalar, Parent=self)


class LPmode(BaseDetector):

    def __init__(self,
                 Mode:          tuple,
                 Source:        Source,
                 Orientation:   str   = 'h',
                 NA:            float = 0.2,
                 Samples:       int   = 401,
                 GammaOffset:   float = 0,
                 PhiOffset:     float = 0,
                 Filter:        float =  'None',
                 Name:          str   = 'Amplitude detector'):

        assert Orientation in ['v','h'], "Orientation should either be v [vertical] or h [horizontal]'"

        self._name, self._coupling = Name, 'Amplitude'

        self._NA = NA

        self.MaxAngle = NA2Angle(self._NA)/180*np.pi

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.ModeNumber = Mode[0]+1, Mode[1]

        self.Source, self.Samples, self.Npts = Source, Samples, int(Samples/8)*2 + 1

        self.Orientation = Orientation

        self.GetFarField()

        self.GetSpherical()



    def GetFarField(self):

        Fiber, CoreDiameter = SMF28()

        temp = fibermodes.field.Field(Fiber.source,
                                      fibermodes.Mode(fibermodes.ModeFamily.HE, *self.ModeNumber),
                                      self.Source.Wavelength,
                                      CoreDiameter*self.Npts/8,
                                      self.Npts).Ex()

        temp = np.array(temp, copy=False)

        temp /= (temp.__abs__()).sum()

        if self.Orientation == 'h': temp = temp.T

        temp = np.fft.fft2(temp)

        temp /= self.GenShift()

        self.Cartesian  = np.fft.fftshift(temp)

        self.Cartesian  /= (self.Cartesian.__abs__()).sum()

        MaxAngle = NA2Angle(self._NA)/180*np.pi

        self.Meshes = AngleMeshes(MaxAngle    = MaxAngle,
                                  Samples     = self.Samples,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)


    def GenShift(self):

        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(self.Npts)*(self.Npts-1)/self.Npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T




    def GetSpherical(self):


        polarImageReal, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.real,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = 40,
                                                                        finalAngle    = 2*np.pi)

        polarImageimag, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.imag,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = 40,
                                                                        finalAngle    = 2*np.pi)

        Scalar = polarImageReal + complex(0,1) * polarImageimag

        self.Scalar = ScalarFarField(Scalar, Parent=self)


        shape = Scalar.imag.shape

        offset = np.pi/2 - self.MaxAngle

        ThetaMesh, PhiMesh = np.mgrid[-np.pi: np.pi:complex(shape[0]),
                                      offset: offset + self.MaxAngle:complex(shape[1])]

        self.Samples = self.Meshes.Phi.Radian.shape

        ZReal = griddata((ThetaMesh.flatten(), PhiMesh.flatten()),
                          Scalar.real.flatten(),
                          (self.Meshes.base[2], self.Meshes.base[1]),
                          fill_value = 0,
                          method     = 'cubic')


        ZImag = griddata((ThetaMesh.flatten(), PhiMesh.flatten()),
                          Scalar.imag.flatten(),
                          (self.Meshes.base[2], self.Meshes.base[1]),
                          fill_value = 0,
                          method     = 'cubic')


        Scalar = ZReal + complex(0,1) * ZImag

        self.Scalar = ScalarFarField(Scalar, Parent=self)































# --
