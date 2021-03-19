#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from PyMieSim.Representations import S1S2, SPF, Stokes, ScalarFarField, Footprint
from PyMieSim.Physics import _Polarization, Angle
from PyMieSim.utils import InterpFull, NA2Angle, GetFieldBinding, Cart2Sp
from PyMieSim.Mesh import FibonacciMesh
from PyMieSim._Coupling import Coupling


class BaseSource(object):

    def __init__(self,
                 Wavelength,
                 Polarization,
                 NA = 0.2):

        self.Wavelength = Wavelength
        self.k = 2 * np.pi / Wavelength
        self.Polarization = _Polarization(Polarization)


class MeshProperty(object):
    """Base class for :class:`Detector` class used to define the properties
    of the angular mesh for Far-Field computations.

    """

    def __init__(self):
        pass

    @property
    def Filter(self):
        return self._Filter

    @Filter.setter
    def Filter(self, val):
        self._Filter = Polarization(val)

    @property
    def PhiOffset(self):
        return self.Mesh.PhiOffset

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.Mesh.UpdateSphere(PhiOffset = val)
        self.GetSpherical()

    @property
    def GammaOffset(self):
        return self.Mesh.GammaOffset

    @GammaOffset.setter
    def GammaOffset(self, val):
        self.Mesh.UpdateSphere(GammaOffset = val)
        self.UnstructuredFarField()

    @property
    def NA(self):
        return self._NA

    @NA.setter
    def NA(self, val):
        if val >= 0.99: val = 0.99
        if val <= 0.01: val = 0.01
        self.MaxAngle = NA2Angle(val).Radian
        self.Mesh.UpdateSphere(MaxAngle = self.MaxAngle)



class BaseDetector(object):
    """Base class for :class:`Detector` class used to define the properties
    of the angular mesh for Far-Field computations.

    """
    def Coupling(self, Scatterer):
        """Return the value of the scattererd light coupling as computed as:

        :math:`|\\iint_{\\Omega}  \Phi_{det} \,\, \\Psi_{scat}^* \,  d \\Omega|^2`

        where :math:`\Phi_{det}` is the capturing field of the detector and
        :math:`\Psi_{scat}` is the scattered field.

        Parameters
        ----------
        Scatterer : :class:`Scatterer`
            Scatterer instance (sphere, cylinder, ...).

        Returns
        -------
        :class:`float`
            Value of the coupling.

        """
        return Coupling(Scatterer = Scatterer, Detector = self)


    def Footprint(self, Scatterer, Num = 200):
        """Return the footprint of the scattererd light coupling with the
        detector as computed as:

        :math:`\\big| \\mathscr{F}^{-1} \\big\\{ \\tilde{ \\psi } (\\xi, \\nu),\
         \\tilde{ \\phi}_{l,m}(\\xi, \\nu)  \\big\\}(\\delta_x, \\delta_y) \\big|^2`.

        where :math:`\\Phi_{det}` is the capturing field of the detector and
        :math:`\\Psi_{scat}` is the scattered field.

        Parameters
        ----------
        Scatterer : :class:`Scatterer`.
            Scatterer instance (sphere, cylinder, ...).

        Returns
        -------
        :class:`Footprint`.
            Dictionnary subclass with all pertienent information.

        """
        return Footprint(Scatterer = Scatterer, Detector = self, Num=Num)


    def SphericalMesh(self,
                      Sampling,
                      MaxAngle,
                      PhiOffset   = 0,
                      GammaOffset = 0,
                      Structured  = True):
        """Method that return an angular mesh (:math:`\\theta`, :math:`\\phi`)
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

        from PyMieSim.utils import PlotUnstructured

        PlotUnstructured(self.Scalar, self.Mesh, Name='Mode field')



class EfficienciesProperties(object):


    @property
    def Qext(self):
        """Extinction efficiency:
        :math:`Q_{ext}=\\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1) / \\text{real} \{ a_n+b_n \}`

        """
        if self._Qext:
            return self._Qext
        else:
            self.GetEfficiencies()
            return self._Qext


    @property
    def Qsca(self):
        """Scattering efficiency:
        :math:`Q_{sca}=\\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)`

        """
        if self._Qsca:
            return self._Qsca
        else:
            self.GetEfficiencies()
            return self._Qsca

    @property
    def Qabs(self):
        """Absorption efficiency:
        :math:`Q_{abs}=Q_{ext}-Q_{sca}`

        """
        if self._Qabs:
            return self._Qabs
        else:
            self.GetEfficiencies()
            return self._Qabs



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
        pass

    def GetEfficiencies(self):
        """Methode compute all Efficiences (:math:`Q_{sca}, Q_{ext}, Q_{abs}`)
        for the scatterer.

        """

        self._Qsca, self._Qext, self._Qabs = self.Bind.Efficiencies



    def S1S2(self, Num=200):
        """Method compute :math:`S_1(\\phi)` and :math:`S_2(\\phi)`.
        For spherical Scatterer such as here S1 and S2 are computed as follow:

        :math:`S_1=\\sum\\limits_{n=1}^{n_{max}} \\frac{2n+1}{n(n+1)}(a_n \\pi_n+b_n \\tau_n)`

        :math:`S_2=\\sum\\limits_{n=1}^{n_{max}}\\frac{2n+1}{n(n+1)}(a_n \\tau_n+b_n \\pi_n)`

        Parameters
        ----------
        Num : :class:`int`
            Number of point (:math:`\\phi`) to evaluate :math:`S_1` and :math:`S_2`.

        Returns
        -------
        :class:`S1S2`
            Description of returned object.

        """

        return S1S2(Parent=self, Num=Num)



    def FarField(self, Num: int = 200):
        """Method Compute scattering Far Field.

        Fields = :math:`\\Big\{ E_{||}(\\phi,\\theta)^2, E_{\\perp}(\\phi,\\theta)^2 \\Big\}`.

        The Fields are up to a constant phase value: :math:`\\exp{(-i k r )}`

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        Returns
        -------
        :class:`ScalarField`
            Dictionnay sub-class with all pertient information as keys.

        """
        return ScalarFarField(Num = Num, Parent = self)



    def _FarField(self, Phi, Theta):
        """Method Compute scattering Far Field for unstructured coordinate.

        Fields = :math:`E_{||}(\\phi,\\theta)^2, E_{\\perp}(\\phi,\\theta)^2`

        The Fields are up to a constant phase value.

        :math:`\\exp{\big(-i k r \big)}`


        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        """

        return self.Bind.UFields(Phi = Phi, Theta=Theta, R=1.)




    def SPF(self, Num=100):
        """Scattering phase function.

        SPF = :math:`E_{\\parallel}(\\phi,\\theta)^2 + E_{\\perp}(\\phi,\\theta)^2`

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\theta , \\phi`) evaluate the SPF [Num, Num].

        Returns
        -------
        :class:`SPF` in :mod:`PyMieSim`
            Dictionnay subclass with all pertient information as keys.

        """

        return SPF(Parent=self, Num=Num)








# -
