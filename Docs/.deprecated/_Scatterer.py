#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import gamma
from PyMieSim.utils import PlotFarField
from PyMieSim.BaseClasses import BaseScatterer, EfficienciesProperties, BaseSource
from PyMieSim.Representations import S1S2, SPF, Field, Stokes
from PyMieSim.Special import Psi, Psi_p, Xi, Xi_p
from PyMieSim.GLMT.sphere import an, bn, cn, dn


class Sphere(BaseScatterer, EfficienciesProperties):
    """Short summary.

    Parameters
    ----------
    Diameter : :class:`float`
        Diameter of the single scatterer in unit of meter.
    Source : :class:`BaseSource`
        Light source object containing info on polarization and wavelength.
    Index : :class:`float`
        Refractive index of scatterer

    Attributes
    ----------
    Area : :class:`float`
        Mathematical 2D area of the scatterer [:math:`\\pi r^2`].
    SizeParam : :class:`float`
        Size parameter of the scatterer [:math:`k r`].

    """

    def __init__(self,
                 Diameter:    float,
                 Source:      BaseSource,
                 Index:       float,
                 IndexMedium: float  = 1.0,
                 MuSphere:    float  = 1.0,
                 MuMedium:    float  = 1.0):

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.nMedium = IndexMedium

        self.Area = np.pi * (Diameter/2)**2

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Qsca, self.Q_ext, self._Qabs = None, None, None

        self.MuSp = MuSphere

        self.Mu = MuMedium

        self._an, self._bn, self._cn, self._dn = [], [], [], []


    def an(self, MaxOrder):
        """ Compute :math:`a_n` coefficient as defined in Eq:III.88 of B&B:

        :math:`a_n = \\frac{
        \mu_{sp} \\Psi_n(\\alpha) \\Psi_n^\prime(\\beta) -
        \\mu M \Psi_n^\prime(\\alpha) \\Psi_n(\\beta)}
        {\mu_{sp} \\xi_n(\\alpha) \\Psi_n^\prime(\\beta)-
        \\mu M \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}`

        With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """

        alpha = self.SizeParam
        beta = self.Index * alpha
        M = self.Index/self.nMedium; MuSp = self.MuSp; Mu = self.Mu;

        for order in range(1, MaxOrder+1):
            numerator = MuSp * Psi(order, alpha) * Psi_p(order, beta) \
                      - Mu * M * Psi_p(order, alpha) * Psi(order, beta)

            denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) \
                        - Mu * M * Xi_p(order, alpha) * Psi(order, beta)

            self._an.append(numerator/denominator)

        return self._an


    def bn(self, MaxOrder):
        """ Compute :math:`b_n` coefficient as defined in Eq:III.89 of B&B:

        :math:`b_n = \\frac{
        \mu M \\Psi_n(\\alpha) \\Psi_n^\prime(\\beta) -
        \\mu_{sp} \Psi_n^\prime(\\alpha) \Psi_n(\\beta)}
        {\mu M \\xi_n(\\alpha) \\Psi_n^\prime(\\beta)-
        \\mu_{sp} \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}`

        With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """
        alpha = self.SizeParam
        beta = self.Index * alpha
        M = self.Index/self.nMedium; MuSp = self.MuSp; Mu = self.Mu;

        for order in range(1, MaxOrder+1):
            numerator = Mu * M * Psi(order, alpha) * Psi_p(order, beta) \
                      - MuSp * Psi_p(order, alpha) * Psi(order, beta)

            denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) \
                        - MuSp  * Xi_p(order, alpha) * Psi(order, beta)

            self._bn.append(numerator/denominator)

        return self._bn


    def cn(self, MaxOrder):
        """ Compute :math:`c_n` coefficient as defined in Eq:III.90 of B&B:

        :math:`c_n = \\frac{
        \mu_{sp} M \\big[ \\xi_n(\\alpha) \\Psi_n^\prime(\\alpha) -
        \\xi_n^\prime(\\alpha) \\Psi_n(\\alpha) \\big]}
        {\mu_{sp} \\xi_n(\\alpha) \\Psi_n^\\prime(\\beta)-
        \\mu M \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}`

        With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """
        alpha = self.SizeParam
        beta = self.Index * alpha
        M = self.Index/self.nMedium; MuSp = self.MuSp; Mu = self.Mu;

        for order in range(1, MaxOrder+1):
            numerator = M * MuSp * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) )

            denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta)

            self._cn.append(numerator/denominator)

        return self._cn


    def dn(self, MaxOrder):
        """ Compute :math:`d_n` coefficient as defined in Eq:III.91 of B&B:

        :math:`d_n = \\frac{
        \mu M^2 \\big[ \\xi_n(\\alpha) \\Psi_n^\prime(\\alpha) -
        \\xi_n^\prime(\\alpha) \\Psi_n(\\alpha) \\big]}
        {\mu M \\xi_n(\\alpha) \\Psi_n^\prime(\\beta)-
        \\mu_{sp} M \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}`

        With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """
        alpha = self.SizeParam
        beta = self.Index * alpha
        M = self.Index/self.nMedium; MuSp = self.MuSp; Mu = self.Mu;

        for order in range(1, MaxOrder+1):

            numerator = Mu * M**2 * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) )

            denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp * Xi_p(order, alpha) * Psi(order, beta)

            self._dn.append(numerator/denominator)

        return self._dn



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


    """
    def __init__(self,
                 g:       float,
                 lc:      float,
                 D:       float,
                 Nc:      float,
                 Source:  BaseSource):

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

        fig0 = PlotFarField(Phi     = Phi,
                            Theta   = Theta,
                            Scalar  = Scalar.reshape([num,num]),
                            Mesh  = scatter,
                            scatter = False,
                            Name    = 'Scattered field')




# -
