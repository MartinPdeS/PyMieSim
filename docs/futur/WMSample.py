#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy       as np
from beartype      import beartype
from scipy.special import gamma
from typing        import Union

from PyMieSim.Tools.units           import Area
from PyMieSim.Source                import PlaneWave, GaussianBeam
from PyMieSim.Tools.Representations import S1S2, SPF, Stokes
from PyMieSim.Tools.ErrorMsg        import *
from PyMieSim.bin.GLMTScatterer     import SPHERE as GLMTSPHERE, CYLINDER as GLMTCYLINDER
from PyMieSim.bin.LMTScatterer      import SPHERE, CYLINDER, SHELLSPHERE1

from PyMieSim.Tools.BaseClasses     import ( BaseScatterer,
                                             ScattererProperties,
                                             BaseSource )


ScalarType = Union[int, float, complex]
SourceType = Union[PlaneWave, GaussianBeam]

class WMSample(object):
    """
    .. note::
        Class representing sample described by the Whittle-Matern RI
        correlation function and using the first Born approximation .

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


    """
    def __init__(self,
                 g      : float,
                 lc     : float,
                 D      : float,
                 Nc     : float,
                 Source : BaseSource):

        self.g  = g; self.lc = lc; self.D  = D; self.Nc = Nc

        self.Source = Source

        self._Perpendicular, self._Parallel = None, None


    def FarField(self, Phi, Theta):

        k = self.Source.k

        term0 = 2 * self.Nc * self.lc * gamma(self.D/2) / np.sqrt(np.pi) * k**4

        term1 = (1-np.sin(Phi-np.pi/2)**2*np.cos(Theta + self.Source.Polarization.Radian)**2)

        term2 = (1 + (2* k * self.lc * np.sin((Phi-np.pi/2)/2)**2)**(self.D/2))

        return term0 * term1 / term2


    def Plot(self, num=200, scatter=False):

        Theta, Phi = np.mgrid[0:2*np.pi:complex(num), -np.pi/2:np.pi/2:complex(num)]

        Scalar = self.GetField(Phi, Theta+np.pi/2)

        from PyMieSim.utils import PlotUnstructured

        PlotUnstructuredAmplitude(Scalar = Scalar, Phi=Phi, Theta=Theta, Name='Mode field')



# -
