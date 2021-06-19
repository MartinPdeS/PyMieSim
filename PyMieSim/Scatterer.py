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
from PyMieSim.GLMT.Scatterer        import ( SPHERE as G_SPHERE,
                                             CYLINDER as G_CYLINDER )

from PyMieSim.LMT.Scatterer   import ( SPHERE,
                                       CYLINDER,
                                       SHELLSPHERE1 )

from PyMieSim.Tools.BaseClasses     import ( BaseScatterer,
                                             ScattererProperties,
                                             BaseSource )


ScalarType = Union[int, float, complex]
SourceType = Union[PlaneWave, GaussianBeam]

class Sphere(BaseScatterer, ScattererProperties):
    """
    .. note::
        Class representing a homogeneous spherical scatterer.

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
        .. note:: Mathematical 2D area of the scatterer [:math:`\\pi r^2`].
    SizeParam : :class:`float`
        .. note:: Size parameter of the scatterer [:math:`k r`].

    """
    kwargformat = [ 'Diameter',
                    'Index',
                    'Material',
                    'nMedium']


    #@beartype
    def __init__(self,
                 Diameter      : float,
                 Source        : SourceType,
                 Index         : ScalarType = None,
                 nMedium       : ScalarType = 1.0,
                 Concentration : ScalarType = None,
                 Material                   = None):

        super().__init__()

        if Material:
            assert Index is None, Error_IndexMaterial
            self.Material    = Material
            self.Index       = Material.Evaluate(Source.Wavelength)

        if Index:
            assert Material is None, Error_IndexMaterial
            self.Material    = None
            self.Index       = Index


        self.type           = 'Sphere'
        self.Diameter       = Diameter
        self.Source         = Source
        self._Concentration = Concentration
        self.nMedium        = nMedium.real#.astype(complex)
        self.Area           = Area(np.pi * (Diameter/2)**2)
        self.SizeParam      = Source.k * ( self.Diameter / 2 )

        self.GetBinding()


    def GetBinding(self):
        """
        .. note::
            Method call and bind c++ scatterer class
        """
        if self.Source.GLMT is True:
            if self.Source._BSC_ is None:
                raise Exception( ErrorGLMTNoBSC )

            self.Bind = G_SPHERE(Index        = self.Index,
                                 Diameter     = self.Diameter,
                                 Wavelength   = self.Source.Wavelength,
                                 nMedium      = self.nMedium,
                                 Polarization = self.Source.Polarization.Radian,
                                 E0           = self.Source.E0,
                                 BSC          = self.Source._BSC_)

        else:
            self.Bind = SPHERE(Index        = self.Index,
                               Diameter     = self.Diameter,
                               Wavelength   = self.Source.Wavelength,
                               nMedium      = self.nMedium,
                               Polarization = self.Source.Polarization.Radian,
                               E0           = self.Source.E0)


    def an(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`a_n` coefficient as defined in Eq:III.88 of B&B:

            .. math::
                a_n = \\frac{
                \mu_{sp} \\Psi_n(\\alpha) \\Psi_n^\prime(\\beta) -
                \\mu M \Psi_n^\prime(\\alpha) \\Psi_n(\\beta)}
                {\mu_{sp} \\xi_n(\\alpha) \\Psi_n^\prime(\\beta)-
                \\mu M \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}

            With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.an(MaxOrder)


    def bn(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`b_n` coefficient as defined in Eq:III.89 of B&B:

            .. math::
                b_n = \\frac{
                \mu M \\Psi_n(\\alpha) \\Psi_n^\prime(\\beta) -
                \\mu_{sp} \Psi_n^\prime(\\alpha) \Psi_n(\\beta)}
                {\mu M \\xi_n(\\alpha) \\Psi_n^\prime(\\beta)-
                \\mu_{sp} \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}

            With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.bn(MaxOrder)


    def cn(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`c_n` coefficient as defined in Eq:III.90 of B&B:

            .. math::
                c_n = \\frac{
                \mu_{sp} M \\big[ \\xi_n(\\alpha) \\Psi_n^\prime(\\alpha) -
                \\xi_n^\prime(\\alpha) \\Psi_n(\\alpha) \\big]}
                {\mu_{sp} \\xi_n(\\alpha) \\Psi_n^\\prime(\\beta)-
                \\mu M \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}

            With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.cn(MaxOrder)


    def dn(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`d_n` coefficient as defined in Eq:III.91 of B&B:

            .. math::
                d_n = \\frac{
                \mu M^2 \\big[ \\xi_n(\\alpha) \\Psi_n^\prime(\\alpha) -
                \\xi_n^\prime(\\alpha) \\Psi_n(\\alpha) \\big]}
                {\mu M \\xi_n(\\alpha) \\Psi_n^\prime(\\beta)-
                \\mu_{sp} M \\xi_n^\\prime (\\alpha) \\Psi_n(\\beta)}

            With :math:`M = \\frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.dn(MaxOrder)



class ShellSphere(BaseScatterer, ScattererProperties):
    """
    .. note::
        Class representing a core+shell spherical scatterer.

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
        .. note:: Mathematical 2D area of the scatterer [:math:`\\pi r^2`].
    SizeParam : :class:`float`
        .. note:: Size parameter of the scatterer [:math:`k r`].

    """

    kwargformat = [ 'CoreDiameter',
                    'ShellWidth',
                    'CoreIndex',
                    'ShellIndex',
                    'nMedium',
                    'Material']
    #@beartype
    def __init__(self,
                 CoreDiameter  : float,
                 ShellWidth    : float,
                 Source        : SourceType,
                 CoreIndex     : ScalarType,
                 ShellIndex    : ScalarType,
                 nMedium       : ScalarType   = 1.0,
                 Concentration : ScalarType   = None,
                 CoreMaterial                 = None,
                 ShellMaterial                = None,
                 ):

        super().__init__()

        ShellDiameter       = CoreDiameter + ShellWidth
        self._Concentration = Concentration

        if all([ CoreIndex, CoreMaterial ] ) or all([ ShellIndex, ShellMaterial ] ) :
            raise AssertionError( Error_IndexMaterial )

        if all([CoreMaterial, ShellMaterial]):
            self.Material = Material
            self.CoreIndex     = CoreMaterial.Evaluate(Source.Wavelength)
            self.ShellIndex    = ShellMaterial.Evaluate(Source.Wavelength)

        if all([CoreIndex, ShellIndex]):
            self.Material = None
            self.CoreIndex     = CoreIndex
            self.ShellIndex    = ShellIndex

        self.type          = '2-Layer Sphere'
        self.CoreDiameter  = CoreDiameter
        self.ShellDiameter = ShellDiameter
        self.Source        = Source
        self.nMedium       = nMedium
        self.Area          = Area(np.pi * (ShellDiameter/2)**2)

        self.GetBinding()


    def GetBinding(self):
        """
        .. note::
            Method call and bind c++ scatterer class
        """
        if self.Source.GLMT is True:
            raise Exception("""This scatterer is not available in the\
                               GLMT framework for now.""")

        else:
            self.Bind = SHELLSPHERE1(ShellIndex     = self.ShellIndex,
                                     CoreIndex      = self.CoreIndex,
                                     ShellDiameter  = self.ShellDiameter,
                                     CoreDiameter   = self.CoreDiameter,
                                     Wavelength     = self.Source.Wavelength,
                                     nMedium        = self.nMedium,
                                     Polarization   = self.Source.Polarization.Radian,
                                     E0             = self.Source.E0)


    def an(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`a_n` coefficient

        """
        return self.Bind.an(MaxOrder)


    def bn(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`b_n` coefficient.

        """
        return self.Bind.bn(MaxOrder)





class Cylinder(BaseScatterer, ScattererProperties):
    """
    .. note::
        Class representing a right angle cylindrical scatterer.

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
        .. note:: Mathematical 2D area of the scatterer [:math:`\\pi r^2`].
    SizeParam : :class:`float`
        .. note:: Size parameter of the scatterer [:math:`k r`].

    """

    kwargformat = [ 'Diameter',
                    'Index',
                    'Material',
                    'nMedium']


    def __init__(self,
                 Diameter      : float,
                 Source        : SourceType,
                 Index         : ScalarType,
                 IndexMedium   : ScalarType  = 1.0,
                 Concentration : ScalarType  = None,
                 ):

        super().__init__()
        self.type           = 'Cylinder'
        self.Diameter       = Diameter
        self.Source         = Source
        self.Index          = Index
        self.nMedium        = IndexMedium
        self.Area           = Area(np.pi * (Diameter/2)**2)
        self.SizeParam      = Source.k * ( self.Diameter / 2 )
        self._Concentration = Concentration

        self.GetBinding()


    def GetBinding(self):
        """
        .. note::
            Method call and bind c++ scatterer class
        """
        if self.Source.GLMT is True:
            if self.Source._BSC_ is None:
                raise Exception( ErrorGLMTNoBSC )

            self.Bind = G_CYLINDER(Index        = self.Index,
                                   Diameter     = self.Diameter,
                                   Wavelength   = self.Source.Wavelength,
                                   nMedium      = self.nMedium,
                                   Polarization = self.Source.Polarization.Radian,
                                   E0           = self.Source.E0,
                                   BSC          = self.Source._BSC_)

        else:
            self.Bind = CYLINDER(Index        = self.Index,
                                 Diameter     = self.Diameter,
                                 Wavelength   = self.Source.Wavelength,
                                 nMedium      = self.nMedium,
                                 Polarization = self.Source.Polarization.Radian,
                                 E0           = self.Source.E0)



    def an(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`a_n` coefficient as defined ref[5]:

            .. math::
                a_n = \\frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
                { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

            | With :math:`m` being the refractive index of the medium and
            |      :math:`m_t` being the refractive index of the index.

        """
        return self.Bind.an(MaxOrder)


    def bn(self, MaxOrder=5):
        """
        .. note::
            Compute :math:`b_n` coefficient as defined in ref[5]:

            .. math::
                b_n = \\frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
                { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

            | With :math:`m` being the refractive index of the medium and
            |      :math:`m_t` being the refractive index of the index.

        """
        return self.Bind.bn(MaxOrder)



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
