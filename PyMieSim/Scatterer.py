import numpy as np
from scipy.special import gamma

from PyMieSim.Mesh import FibonacciMesh
from PyMieSim.utils import PlotFarField
from PyMieSim.Physics import Angle, Source
from PyMieSim.BaseClasses import BaseScatterer, EfficienciesProperties
from PyMieSim.Representations import S1S2, SPF, Field, Stokes



class Sphere(BaseScatterer, EfficienciesProperties):
    """Short summary.

    Parameters
    ----------
    Diameter : :class:`float`
        Diameter of the single scatterer in unit of meter.
    Source : :class:`Source`
        Light source object containing info on polarization and wavelength.
    Index : :class:`float`
        Refractive index of scatterer

    Attributes
    ----------
    Area : :class:`float`
        Mathematical 2D area of the scatterer [:math:`\\pi r^2`].
    SizeParam : :class:`float`
        Size parameter of the scatterer [:math:`k r`].
    _Stokes : :class:`Stokes`
        Stoke representation class
    _SPF : :class:`SPF`
        Scattering phase function representation class
    _Parallel : :class:`Field`
        Parallel field representation class
    _Perpendicular : :class:`Field`
        Perpendicular field representation class
    _S1S2 : :class:`S1S2`
        S1 and S2 values representation class
    _phi : :class:`list`
        Last phi list used for computing S1S2 or Field or SPF or Stokes
    _theta : :class:`list`
        Last theta list used for computing S1S2 or Field or SPF or Stokes

    """

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 IndexMedium: float  = 1.0):

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.nMedium = IndexMedium

        self.Area = np.pi * (Diameter/2)**2

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        self._phi, self._theta = [None], [None]

        self._Qsca, self.Q_ext, self._Qabs = None, None, None




class WMSample(object):
    """Sample represented by the Whittle-Matern RI correlation function and
    using the first Born approximation .

    Parameters
    ----------
    g : :class:`float`
        Description of parameter `g`.
    lc : :class:`float`
        Correlation lenght of RI of the sample
    D : :class:`float`
        Form factor of the sample
    Nc : :class:`float`
        Scalling factor of the sample.
    Source : :class:`Source`
        Light source object containing info on polarization and wavelength.
    Mesh : :class:`FibonacciMesh`
        Object FibonacciMesh which describes the point in fourier space that are
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



    def Parallel(self, Mesh):
        if not isinstance(self._Parallel, np.ndarray):
            self.GetField(Mesh.Phi.Radian, Mesh.Theta.Radian)
            return self._Parallel
        else:
            return self._Parallel


    def Perpendicular(self, Mesh):
        if not isinstance(self._Perpendicular, np.ndarray):
            self.GetField(Mesh.Phi.Radian, Mesh.Theta.Radian)
            return self._Perpendicular*0
        else:
            return self._Perpendicular


    def GetField(self, Phi, Theta):

        k = self.Source.k

        term0 = 2 * self.Nc * self.lc * gamma(self.D/2) / np.sqrt(np.pi) * k**4

        term1 = (1-np.sin(Phi-np.pi/2)**2*np.cos(Theta + self.Source.Polarization.Radian)**2)

        term2 = (1 + (2* k * self.lc * np.sin((Phi-np.pi/2)/2)**2)**(self.D/2))

        return term0 * term1 / term2


    def Plot(self, num=200, scatter=False):

        Theta, Phi = np.mgrid[0:2*np.pi:complex(num), -np.pi/2:np.pi/2:complex(num)]

        Scalar = self.GetField(Phi, Theta+np.pi/2)

        fig0 = PlotFarField(Phi     = Phi,
                            Theta   = Theta,
                            Scalar  = Scalar.reshape([num,num]),
                            Mesh  = scatter,
                            scatter = False,
                            Name    = 'Scattered field')








# -
