import numpy as np
from scipy.special import gamma

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.utils import Source, PlotFarField
from PyMieCoupling.classes.BaseClasses import BaseScatterer





class Scatterer(BaseScatterer):
    """Single scatterer object.

    Parameters
    ----------
    Diameter : float
        Diameter of the single scatterer in unit of meter.
    Source : Source
        Light source object containing info on polarization and wavelength.
    Index : float
        Refractive index of scatterer
    Meshes : AngleMeshes
        Object AngleMeshes which describes the point in fourier space that are
        used for computation.

    """

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,) -> None:

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        self._phi, self._theta = [None], [None]







class FullScatterer(BaseScatterer):
    """Full-field Single scatterer object.

    Parameters
    ----------
    Diameter : float
        Diameter of the single scatterer in unit of meter.
    Source : Source
        Light source object containing info on polarization and wavelength.
    Index : float
        Refractive index of scatterer

    """
    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Sampling:    int     = 1000):

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        self.Meshes = AngleMeshes(MaxAngle    = np.deg2rad(180),
                                  Sampling    = Sampling,
                                  PhiOffset   = 0,
                                  GammaOffset = 0)





class WMSample(object):
    """Sample represented by the Whittle-Matern RI correlation function and
    using the first Born approximation .

    Parameters
    ----------
    g : float
        Description of parameter `g`.
    lc : float
        Correlation lenght of RI of the sample
    D : float
        Form factor of the sample
    Nc : float
        Scalling factor of the sample.
    Source : Source
        Light source object containing info on polarization and wavelength.
    Meshes : AngleMeshes
        Object AngleMeshes which describes the point in fourier space that are
        used for computation.

    """
    def __init__(self,
                 g:       float,
                 lc:      float,
                 D:       float,
                 Nc:      float,
                 Source:  Source):

        self.g  = g; self.lc = lc; self.D  = D; self.Nc = Nc

        self.Source = Source

        self._Perpendicular, self._Parallel = None, None



    def Parallel(self, Meshes):
        if not isinstance(self._Parallel, np.ndarray):
            self.GetScalar(Meshes)
            return self._Parallel
        else:
            return self._Parallel

    #@property
    def Perpendicular(self, Meshes):
        if not isinstance(self._Perpendicular, np.ndarray):
            self.GetScalar(Meshes)
            return self._Perpendicular*0
        else:
            return self._Perpendicular


    def GetField(self, Phi, Theta):

        k = self.Source.k

        term0 = 2 * self.Nc * self.lc * gamma(self.D/2) / np.sqrt(np.pi) * k**4

        term1 = (1-np.sin(Phi-np.pi/2)**2*np.cos(Theta + self.Source.Polarization.Radian)**2)

        term2 = (1 + (2* k * self.lc * np.sin((Phi-np.pi/2)/2)**2)**(self.D/2))

        return term0 * term1 / term2


    def GetScalar(self, Meshes):

        return self.GetField(Meshes.Phi.Radian, Meshes.Theta.Radian)


    def Plot(self, num=200, scatter=False):

        Theta, Phi = np.mgrid[0:2*np.pi:complex(num), -np.pi/2:np.pi/2:complex(num)]

        Scalar = self.GetField(Phi, Theta+np.pi/2)

        fig0 = PlotFarField(Phi     = Phi,
                            Theta   = Theta,
                            Scalar  = Scalar.reshape([num,num]),
                            Meshes  = scatter,
                            scatter = False,
                            Name    = 'Scattered field')





# -
