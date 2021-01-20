
import numpy as np
import fibermodes

from PyMieCoupling.classes.BaseClasses import BaseDetector
from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Source, SMF28, Angle, _Polarization, PlotUnstructureData, interp_at, InterpFull
from PyMieCoupling.physics import FraunhoferDiffraction



class Photodiode(BaseDetector):

    def __init__(self,
                 Source:            Source = None,
                 NA:                float  = 0.2,
                 Sampling:           int    = 401,
                 GammaOffset:       float  = 0,
                 PhiOffset:         float  = 0,
                 Filter:            float  = 'None',
                 Name:              str    = 'Intensity Detector'):

        self._name, self._coupling  = Name, 'Intensity'

        self._NA = NA

        self.Source, self.Sampling = Source, Sampling

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.Meshes = AngleMeshes(MaxAngle    = NA2Angle(self._NA).Radian,
                                  Sampling    = self.Sampling,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)

        self.Sampling = self.Meshes.Phi.Radian.shape

        self.GetSpherical()


    def GetSpherical(self):

        Scalar = np.ones(self.Sampling) / (self.Sampling)

        self.Scalar = Scalar


class LPmode(BaseDetector):

    def __init__(self,
                 Mode:           tuple,
                 Source:         Source,
                 Orientation:    str   = 'h',
                 NA:             float = 0.2,
                 Sampling:       int   = 401,
                 InterpSampling: int   = 101,
                 GammaOffset:    float = 0,
                 PhiOffset:      float = 0,
                 Filter:         float =  'None',
                 Name:           str   = 'Amplitude detector'):

        assert Orientation in ['v','h'], "Orientation should either be v [vertical] or h [horizontal]'"

        self._name, self._coupling = Name, 'Amplitude'

        self._NA = NA

        self.MaxAngle = NA2Angle(self._NA).Radian

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter); self.Orientation = Orientation

        self.ModeNumber = Mode[0]+1, Mode[1]

        self.Source, self.Sampling, self.InterpSampling = Source, Sampling, InterpSampling

        self.debug = False

        MaxAngle = NA2Angle(self._NA).Radian

        self.Meshes = AngleMeshes(MaxAngle    = MaxAngle,
                                  Sampling    = self.Sampling,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)

        self.Sampling = self.Meshes.Phi.Radian.shape

        self.GetFarField()

        self.GetSpherical()



    def GetFarField(self):

        Fiber, CoreDiameter = SMF28()

        temp = fibermodes.field.Field(Fiber.source,
                                      fibermodes.Mode(fibermodes.ModeFamily.HE, *self.ModeNumber),
                                      self.Source.Wavelength,
                                      CoreDiameter*self.InterpSampling/4,
                                      self.InterpSampling).Ex()

        temp = np.array(temp, copy=False)

        if self.Orientation == 'h': temp = temp.T

        self.Cartesian = FraunhoferDiffraction(temp)




    def GetSpherical(self):

        shape = self.Cartesian.shape

        ThetaMesh, PhiMesh = np.mgrid[-np.pi: np.pi:complex(shape[0]),
                                      np.pi/2: np.pi/2-self.MaxAngle*1.01 :complex(shape[1])]

        Scalar = interp_at(PhiMesh.flatten(),
                           ThetaMesh.flatten(),
                           self.Cartesian.astype(np.complex).flatten(),
                           self.Meshes.base.Phi,
                           self.Meshes.base.Theta,
                           algorithm='linear',
                           extrapolate=True)

        norm = np.sqrt( np.sum((self.Meshes.SinMesh * Scalar.__abs__())**2) )

        Scalar /=  norm

        if self.debug:
            PlotUnstructureData(Scalar, self.Meshes.base.Theta, self.Meshes.base.Phi)

        self.Scalar = Scalar
