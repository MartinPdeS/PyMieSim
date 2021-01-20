import numpy as np
from scipy.special import gamma

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.utils import Source, PlotFarField
from PyMieCoupling.classes.BaseClasses import BaseScatterer





class Scatterer(BaseScatterer):

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Meshes:      AngleMeshes  = None,
                 Acceptance:  list         = 20,
                 Sampling:     int          = 1000,
                 GammaOffset: float        = 0,
                 PhiOffset:   float        = 0) -> None:

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.Acceptance = np.deg2rad(Acceptance)

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        if Meshes:
            self.Meshes = Meshes
        else:
            self.Meshes = AngleMeshes(MaxAngle    = self.Acceptance,
                                      Sampling     = Sampling,
                                      PhiOffset   = PhiOffset,
                                      GammaOffset = GammaOffset)





class FullScatterer(BaseScatterer):

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Sampling:    int     = 1000):

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.Acceptance = np.deg2rad(180)

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        self.Meshes = AngleMeshes(MaxAngle    = self.Acceptance,
                                  Sampling     = Sampling,
                                  PhiOffset   = 0,
                                  GammaOffset = 0)





class Sample(object):
    def __init__(self,
                 g,
                 lc,
                 D,
                 Nc,
                 Source,
                 Meshes      = None,
                 Acceptance  = 20,
                 Sampling    = 1000,
                 GammaOffset = 0,
                 PhiOffset   = 0) -> None:

        self.g  = g; self.lc = lc; self.D  = D; self.Nc = Nc

        self.Acceptance = Acceptance; self.Source = Source

        if Meshes:
            self.Meshes = Meshes
        else:
            self.Meshes = AngleMeshes(MaxAngle    = self.Acceptance,
                                      Sampling    = Sampling,
                                      PhiOffset   = PhiOffset,
                                      GammaOffset = GammaOffset)

        self.Parallel = self.GetScalar(self.Meshes.Theta.Radian, self.Meshes.Phi.Radian)

        self.Perpendicular = self.Parallel


    def GetScalar(self, Phi, Theta):

        k = self.Source.k

        term0 = 2 * self.Nc * self.lc * gamma(self.D/2) / np.sqrt(np.pi) * k**4

        term1 = (1-np.sin(Phi-np.pi/2)**2*np.cos(Theta + self.Source.Polarization.Radian)**2)

        term2 = (1 + (2* k * self.lc * np.sin((Phi-np.pi/2)/2)**2)**(self.D/2))

        return term0 * term1 / term2


    def Plot(self, num=200, scatter=True):

        Theta, Phi = np.mgrid[0:2*np.pi:complex(num), -np.pi/2:np.pi/2:complex(num)]

        Scalar = self.GetScalar(Phi, Theta+np.pi/2)

        fig0 = PlotFarField(Phi     = Phi,
                            Theta   = Theta,
                            Scalar  = Scalar.reshape([num,num]),
                            Meshes  = self.Meshes,
                            Name    = 'Scattered field',
                            scatter = scatter)





# -
