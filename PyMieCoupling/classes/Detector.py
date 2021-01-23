
import numpy as np
import fibermodes
from ai import cs

from PyMieCoupling.classes.BaseClasses import BaseDetector, MeshProperty
from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Source, SMF28, Angle, _Polarization, PlotUnstructureData, interp_at
from PyMieCoupling.physics import FraunhoferDiffraction



class Photodiode(BaseDetector, MeshProperty):
    """Short summary.

    Parameters
    ----------
    Source : Source
        Light source object containing info on polarization and wavelength.
    NA : float
        Numerical aperture of imaging system.
    Sampling : int
        Number of sampling points for the mode (inside NA).
    GammaOffset : float
        Angle offset of detector in the direction perpendicular to polarization.
    PhiOffset : float
        Angle offset of detector in the direction parallel to polarization.
    Filter : float
        Angle of polarization filter in front of detector. Default is "None"
    CouplingMode : str
        Methode for computing mode coupling. Either Centered or Mean.
    """
    def __init__(self,
                 Source:            Source = None,
                 NA:                float  = 0.2,
                 Sampling:          int    = 401,
                 GammaOffset:       float  = 0,
                 PhiOffset:         float  = 0,
                 Filter:            float  = 'None',
                 CouplingMode:      str    = 'Centered'):


        self._CouplingMode = ('Intensity', CouplingMode)

        self.Source = Source

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.Meshes = AngleMeshes(MaxAngle    = NA2Angle(NA).Radian,
                                  Sampling    = Sampling,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)

        self.GetSpherical()


    def GetSpherical(self):

        self.Scalar = np.ones(self.Meshes.Sampling) / (self.Meshes.Sampling)



class LPmode(BaseDetector, MeshProperty):
    """Short summary.

    Parameters
    ----------
    Mode : tuple
        LP mode index l, m.
    Source : Source
        Light source object containing info on polarization and wavelength.
    NA : float
        Numerical aperture of imaging system.
    Sampling : int
        Number of sampling points for the mode (inside NA).
    InterpSampling : int
        Number of sampling point for interpolation of FarField mode.
    GammaOffset : float
        Angle offset of detector in the direction perpendicular to polarization.
    PhiOffset : float
        Angle offset of detector in the direction parallel to polarization.
    Filter : float
        Angle of polarization filter in front of detector. Default is "None"
    CouplingMode : str
        Methode for computing mode coupling. Either Centered or Mean.
    """

    def __init__(self,
                 Mode:           tuple,
                 Source:         Source,
                 NA:             float = 0.2,
                 Sampling:       int   = 401,
                 InterpSampling: int   = 101,
                 GammaOffset:    float = 0,
                 PhiOffset:      float = 0,
                 Filter:         float =  'None',
                 CouplingMode:   str   = 'Centered',):

        if len(Mode) <= 2: Mode = Mode[0], Mode[1], 'h'

        assert Mode[2] in ['v','h',''], "Mode orientation should either be v [vertical] or h [horizontal]"

        assert CouplingMode in ['Centered','Mean'], "Coupling mode can either be Centered or Mean"

        assert NA < 0.5, "Numerical aperture has to be under 0.5 radian"

        self._CouplingMode = ('Amplitude', CouplingMode)

        self._Filter = _Polarization(Filter)

        self.ModeNumber = Mode[0]+1, Mode[1], Mode[2]

        self.Source, self.InterpSampling = Source, InterpSampling

        self.Meshes = AngleMeshes(MaxAngle    = NA2Angle(NA).Radian,
                                  Sampling    = Sampling,
                                  PhiOffset   = PhiOffset,
                                  GammaOffset = GammaOffset)

        self.GetFarField()

        self.GetSpherical()



    def GetFarField(self):

        Fiber, CoreDiameter = SMF28()

        temp = fibermodes.field.Field(Fiber.source,
                                      fibermodes.Mode(fibermodes.ModeFamily.HE, *self.ModeNumber[:2]),
                                      self.Source.Wavelength,
                                      CoreDiameter*self.InterpSampling/8,
                                      self.InterpSampling).Ex()

        temp = np.array(temp, copy=False)

        if self.ModeNumber == 'h': temp = temp.T

        self.Cartesian = FraunhoferDiffraction(temp)



    def GetSpherical(self):

        shape = self.Cartesian.shape

        x, y = np.mgrid[-50: 50: complex(shape[0]), -50: 50: complex(shape[1])]

        z = 50 / np.tan(self.Meshes.MaxAngle*2)

        r, phi, theta = cs.cart2sp(x.flatten(), y.flatten(), x.flatten()*0+z)

        self.Scalar = interp_at(x           = phi.flatten(),
                                y           = theta.flatten(),
                                v           = self.Cartesian.astype(np.complex).flatten(),
                                xp          = self.Meshes.base.Phi,
                                yp          = self.Meshes.base.Theta,
                                algorithm   = 'linear',
                                extrapolate = True)

        norm = np.sqrt( np.sum((self.Meshes.SinMesh * self.Scalar.__abs__())**2) )

        self.Scalar /=  norm

        if False:
            PlotUnstructureData(Scalar, self.Meshes.base.Theta, self.Meshes.base.Phi)

















# -
