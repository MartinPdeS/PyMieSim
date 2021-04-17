#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from PyMieSim.Representations import S1S2, SPF, Stokes, ScalarFarField, Footprint
from PyMieSim.Physics import _Polarization, Angle
from PyMieSim.utils import InterpFull, NA2Angle, Cart2Sp
from PyMieSim.units import Power, Area
from PyMieSim.Mesh import FibonacciMesh
from PyMieSim._Coupling import Coupling
import PyMieSim.Plots as plot
from PyMieSim.Constants import *

EPS = 1e-6

class BaseSource(object):

    def __init__(self,
                 Wavelength,
                 Polarization,
                 NA = 0.2):

        pass

    @property
    def Polarization(self):
        return self._Polarization

    @Polarization.setter
    def Polarization(self, val):
        self._Polarization = _Polarization(val)


class MeshProperty(object):
    """
    Base class for :class:`Detector` class used to define the properties
    of the angular mesh for Far-Field computations.

    """

    @property
    def Filter(self):
        return self._Filter

    @Filter.setter
    def Filter(self, val):
        self._Filter = _Polarization(val)

    @property
    def PhiOffset(self):
        return self.Mesh.PhiOffset

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.Mesh.UpdateSphere(PhiOffset = val)
        self.FarField(Structured=False)

    @property
    def GammaOffset(self):
        return self.Mesh.GammaOffset

    @GammaOffset.setter
    def GammaOffset(self, val):
        self.Mesh.UpdateSphere(GammaOffset = val)
        self.FarField(Structured=False)

    @property
    def NA(self):
        return self._NA

    @NA.setter
    def NA(self, val):
        if val <= EPS: val = EPS

        self.MaxAngle = NA2Angle(val).Radian
        self.Mesh.UpdateSphere(MaxAngle = self.MaxAngle)


class BaseDetector(object):
    """
    Base class for :class:`Detector` class used to define the properties
    of the angular mesh for Far-Field computations.

    """

    def _Coupling(self, Scatterer):
        """
        Return the value of the scattererd light coupling as computed as:

        .. math::
            |\\iint_{\\Omega}  \Phi_{det} \,\, \\Psi_{scat}^* \,  d \\Omega|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        Parameters
        ----------
        Scatterer : :class:`Scatterer`
            Scatterer instance (sphere, cylinder, ...).

        Returns
        -------
        :class:`float`
            Value of the coupling.

        """

        return Coupling(Scatterer = Scatterer, Detector = self)# * Scatterer.nMedium * c * eps0 / 2


    def Coupling(self, Scatterer):
        """
        Return coupling power which is detected.

        Parameters
        ----------
        Scatterer : :class:`Scatterer`
            Scatterer instance (sphere, cylinder, ...).

        Returns
        -------
        :class:`float`
            Value of the coupling power [:math:`W`].

        """
        C = self._Coupling(Scatterer) * eps0 * c * 0.5

        return Power(C)


    def Footprint(self, Scatterer, Num = 200):
        """
        Return the footprint of the scattererd light coupling with the
        detector as computed as:

        .. math::
            \\big| \\mathscr{F}^{-1} \\big\\{ \\tilde{ \\psi } (\\xi, \\nu),\
                   \\tilde{ \\phi}_{l,m}(\\xi, \\nu)  \\big\\}(\\delta_x, \\delta_y) \\big|^2

        | Where:
        |   :math:`\\Phi_{det}` is the capturing field of the detector and
        |   :math:`\\Psi_{scat}` is the scattered field.

        Parameters
        ----------
        Scatterer : :class:`Scatterer`.
            Scatterer instance (sphere, cylinder, ...).

        Returns
        -------
        :class:`Footprint`.
            Dictionnary subclass with all pertienent information.

        """
        return Footprint(Scatterer = Scatterer, Detector = self, Num = Num)


    def SphericalMesh(self,
                      Sampling,
                      MaxAngle,
                      PhiOffset   = 0,
                      GammaOffset = 0,
                      Structured  = True):
        """
        Method that return an angular mesh (:math:`\\theta`, :math:`\\phi`)
        which is either structured or not. If not the pattern follow a
        Fibonacci mesh.

        """

        if Structured:
            x, y = np.mgrid[-50: 50: complex(Sampling), -50: 50: complex(Sampling)]
            z = 50 / np.tan(MaxAngle)
            _, theta, phi = Cart2Sp(x, y, x*0+z)

            return phi, theta

        else:
            return FibonacciMesh(MaxAngle    = MaxAngle,
                                 Sampling    = Sampling,
                                 PhiOffset   = PhiOffset,
                                 GammaOffset = GammaOffset)


    def Plot(self):
        """Method that plot the real part of the scattered field
        (:math:`E_{\\theta}` and :math:`E_{\\phi}`).

        """

        plot.Unstructured(Scalar = self.Scalar,
                          Mesh   = self.Mesh,
                          Name   = 'Mode field',
                          Mode   = 'Amplitude')


class EfficienciesProperties(object):


    @property
    def Qext(self):
        """
        Extinction efficiency.

        .. math::
            Q_{ext}=\\frac{2}{x^2} \sum\limits_{n=1}^{n_{max}}  (2n+1) / \\text{real} \{ a_n+b_n \}

        """
        if self._Qext:
            return self._Qext
        else:
            self.GetEfficiencies()
            return self._Qext


    @property
    def Qsca(self):
        """
        Scattering efficiency.

        .. math::
            Q_{sca}=\\frac{2}{x^2} \sum\limits_{n=1}^{n_{max}} (2n+1)(|a_n|^2+|b_n|^2)

        """
        if self._Qsca:
            return self._Qsca
        else:
            self.GetEfficiencies()
            return self._Qsca


    @property
    def Qabs(self):
        """
        Absorption efficiency.

        .. math::
            Q_{abs} = Q_{ext}-Q_{sca}

        """
        if self._Qabs:
            return self._Qabs
        else:
            self.GetEfficiencies()
            return self._Qabs


    @property
    def Qback(self):
        """
        Backscattering efficiency.

        .. math::
            Q_{back} = \\frac{1}{x^2} \\Big| \sum\limits_{n=1}^{n_{max}} (2n+1)(-1)^n (a_n - b_n) \\Big|^2

        """
        if self._Qback:
            return self._Qback
        else:
            self.GetEfficiencies()
            return self._Qback


    @property
    def Qratio(self):
        """
        Ratio of backscattering over total scattering.

        .. math::
            Q_{ratio} = \\frac{Q_{back}}{Q_{sca}}

        """
        if self._Qratio:
            return self._Qratio
        else:
            self.GetEfficiencies()
            return self._Qratio

    @property
    def g(self):
        """
        Ratio of backscattering over total scattering.

        .. math::
            g = \\frac{4}{Q_{sca} x^2}
            \\left[ \\sum\\limits_{n=1}^{n_{max}} \\frac{n(n+2)}{n+1} \\text{Re} \\left
            \{a_n a_{n+1}^* + b_n b_{n+1}^*\\right\} + \\sum\\limits_{n=1}^{n_{max}}
            \\frac{2n+1}{n(n+1)} \\text{Re} \\left\{ a_n b_n^* \\right\} \\right]}

        """
        if self._g:
            return self._g
        else:
            self.GetEfficiencies()
            return self._g


    @property
    def Qpr(self):
        """
        Ratio of backscattering over total scattering.

        .. math::
            Q_{pr} = Q_{ext} - g * Q_{sca}

        """
        if self._Qpr:
            return self._Qpr
        else:
            self.GetEfficiencies()
            return self._Qpr


    @property
    def Efficiencies(self):
        """Methode compute all Efficiences (:math:`Q_{sca}, Q_{ext}, Q_{abs},
        Q_{back}, Q_{ratio}, g, Q_{pr}`) for the scatterer.

        """
        from PyMieSim.Representations import Efficiences

        if self._Efficiencies:
            return self._Efficiencies
        else:
            self._Efficiencies = Efficiences(self)
            return self._Efficiencies


class BaseScatterer(object):
    """Base class for :class:`Sphere`.
    This class containes all the methodes that output something interesting for
    the user.

    Parameters
    ----------
    diameter : :class:`float`
        Diameter of the scatterer.
    wavelength : :class:`float`
        Wavelength of the incident lightfield.
    index : :class:`float`
        Refractive index of the scatterer.
    npts : :class:`int`
        Number of points for the full solid angle of the far-field, later to
        be interpolated.


    """

    def __init__(self):
        self._Qsca   = None
        self._Qext   = None
        self._Qabs   = None
        self._Qback  = None
        self._Qratio = None
        self._g      = None
        self._Qpr    = None
        self._an     = []
        self._bn     = []
        self._cn     = []
        self._dn     = []
        self._Efficiencies = None


    def GetEfficiencies(self):
        """Methode compute all Efficiences (:math:`Q_{sca}, Q_{ext}, Q_{abs},
        Q_{back}, Q_{ratio}, g, Q_{pr}`) for the scatterer.

        """

        (self._Qsca,
         self._Qext,
         self._Qabs,
         self._Qback,
         self._Qratio,
         self._g,
         self._Qpr )   = self.Bind.Efficiencies


    def S1S2(self, Num=200):
        """
        Method compute :math:`S_1(\\phi)` and :math:`S_2(\\phi)`.
        For spherical Scatterer such as here S1 and S2 are computed as follow:

        .. math::
            S_1=\\sum\\limits_{n=1}^{n_{max}} \\frac{2n+1}{n(n+1)}(a_n \\pi_n+b_n \\tau_n)

            S_2=\\sum\\limits_{n=1}^{n_{max}}\\frac{2n+1}{n(n+1)}(a_n \\tau_n+b_n \\pi_n)

        Parameters
        ----------
        Num : :class:`int`
            Number of point (:math:`\\phi`) to evaluate :math:`S_1` and :math:`S_2`.

        Returns
        -------
        :class:`dict`
            Dictionnay sub-class with all pertient information as keys.

        """

        return S1S2(Parent=self, Num=Num)


    def Stokes(self, Num=200):
        """
        Method compute and return the Stokes parameters: I, Q, U, V.
        Those parameters are defined as:

        .. math:
            I &= \\big| E_x \big|^2 + \\big| E_y \\big|^2

            Q &= \\big| E_x \big|^2 - \\big| E_y \\big|^2

            U &= 2 \\mathcal{Re} \\big\{ E_x E_y^* \\big\}

            V &= 2 \\mathcal{Im} \\big\{ E_x E_y^* \\big\}

        Parameters
        ----------
        Num : :class:`int`
            Number of point (:math:`\\phi`) to evaluate :math:`S_1` and :math:`S_2`.

        Returns
        -------
        :class:`dict`
            Dictionnay sub-class with all pertient information as keys.

        """

        return Stokes(Parent=self, Num=Num)



    def FarField(self, Num: int = 200):
        """Method Compute scattering Far Field.

        .. math::
            \\text{Fields} = E_{||}(\\phi,\\theta)^2,
                             E_{\\perp}(\\phi,\\theta)^2


        The Fields are up to a constant phase value:

        .. math::
            \\exp{\\big(-i k r \\big)}

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        Returns
        -------
        :class:`dict`
            Dictionnay sub-class with all pertient information as keys.

        """
        return ScalarFarField(Num = Num, Parent = self)



    def uFarField(self, Phi, Theta, R):
        """Method Compute scattering Far Field for unstructured coordinate.

        .. math::
            \\text{Fields} = E_{||}(\\phi,\\theta)^2,
                             E_{\\perp}(\\phi,\\theta)^2


        The Fields are up to a constant phase value.

        .. math::
            \\exp{\\big(-i k r \\big)}


        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        Returns
        -------
        :class:`np.array`
            The unstructured far-field

        """

        return self.Bind.uFields(Phi = Phi, Theta=Theta, R=R)


    def sFarField(self, Phi, Theta, R):
        """Method Compute scattering Far Field for structured coordinate.

        .. math::
            \\text{Fields} = E_{||}(\\phi,\\theta)^2,
                             E_{\\perp}(\\phi,\\theta)^2


        The Fields are up to a constant phase value.

        .. math::
            \\exp{\\big(-i k r \\big)}


        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        Returns
        -------
        :class:`np.array`
            The structured far-field

        """

        return self.Bind.sFields(Phi = Phi, Theta=Theta, R=R)


    def uS1S2(self, Phi, Theta):
        """Method Compute scattering Far Field for unstructured coordinate.

        .. math::
            \\text{Fields} = E_{||}(\\phi,\\theta)^2,
                             E_{\\perp}(\\phi,\\theta)^2

        The Fields are up to a constant phase value.

        .. math::
            \\exp{\\big(-i k r \\big)}


        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        Returns
        -------
        :class:`np.array`
            The unstructured non-propagated far-field


        """

        return self.Bind.uS1S2(Phi = Phi, Theta=Theta)


    def sS1S2(self, Phi, Theta):
        """Method Compute scattering Far Field for structured coordinate.

        .. math::
            \\text{Fields} = E_{||}(\\phi,\\theta)^2,
                             E_{\\perp}(\\phi,\\theta)^2

        The Fields are up to a constant phase value.

        :math:`\\exp{\big(-i k r \big)}`


        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        Returns
        -------
        :class:`np.array`
            The structured non-propagated far-field

        """

        return self.Bind.sS1S2(Phi = Phi, Theta=Theta)


    def SPF(self, Num=100):
        """Scattering phase function.

        .. math::
            \\text{SPF} = E_{\\parallel}(\\phi,\\theta)^2 + E_{\\perp}(\\phi,\\theta)^2

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\theta , \\phi`) evaluate the SPF [Num, Num].

        Returns
        -------
        :class:`dict`
            Dictionnay subclass with all pertinent information as keys.

        """

        return SPF(Parent=self, Num=Num)


    def PoyntingVector(self, Mesh):
        """
        Method return the Poynting vector norm defined as:

        .. math::
            \\vec{S} = \\epsilon c^2 \\vec{E} \\times \\vec{B}

        Parameters
        ----------
        Mesh : :class:`FibonacciMesh`
            Number of voxel in the 4 pi space to compute energy flow.

        Returns
        -------
        :class:`np.array`
            Poynting field [:math:`W/m^2`]

        """

        EPhi, ETheta = self.uFarField(Mesh.Phi.Radian, Mesh.Theta.Radian,1.)

        NormE        = np.sqrt(np.abs(EPhi)**2 + np.abs(ETheta)**2)

        NormB        = NormE/c

        Poynting     = eps0 * c**2 * NormE * NormB    #TODO change eps0 for eps

        return Poynting


    def EnergyFlow(self, Mesh):
        """
        Method return energy flow defined as:

        .. math::

            W_a &= \\sigma_{sca} * I_{inc}

            P &= \\int_{A} I dA

            I &= \\frac{c n \\epsilon_0}{2} |E|^2

        | With:
        |     I : Energy density
        |     n  : Refractive index of the medium
        |     :math:`\\epsilon_0` : Vaccum permitivity
        |     E  : Electric field
        |    :math:`\\sigma_{sca}`: Scattering cross section.

        More info on wikipedia link (see ref[6]).


        Parameters
        ----------
        Mesh : :class:`FibonacciMesh`
            Number of voxel in the 4 pi space to compute energy flow.

        Returns
        -------
        :class:`float`
            Energy flow [:math:`W`]

        """

        Poynting = self.PoyntingVector(Mesh)

        if Mesh.Structured:
            Wtotal = 0.5 * np.sum( Poynting * Mesh.SinMesh ) * Mesh.dOmega.Radian

        else:
            Wtotal = 0.5 * np.sum( Poynting ) * Mesh.dOmega.Radian

        return Power(Wtotal)


    def CrossSection(self, Mesh):
        """
        Method return scattering cross section, see :func:`EnergyFlow`

        Parameters
        ----------
        Mesh : :class:`FibonacciMesh`
            Number of voxel in the 4 pi space to compute scattering cross section.

        Returns
        -------
        :class:`float`
            scattering cross section [:math:`m^2`]
        """
        #return self.EnergyFlow(Mesh) / self.Source.I

        return Area(self.Qsca*self.Area)
















# -
