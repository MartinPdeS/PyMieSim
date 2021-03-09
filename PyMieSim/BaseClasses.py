#!/usr/bin/env python
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import numpy as np
import matplotlib
import cartopy.crs as ccrs
from ai import cs

from PyMieSim.Representations import S1S2, SPF, Stokes, Field, ScalarFarField
from PyMieSim.Representations import Footprint
from PyMieSim.Physics import _Polarization, Angle
from PyMieSim.utils import InterpFull, NA2Angle
from PyMieSim.Mesh import FibonacciMesh
from PyMieSim._Coupling import Coupling
from PyMieSim.LMT.Sphere import S1S2 as LMTS1S2, Fields as LMTFields
from PyMieSim.GLMT.Sphere import S1S2 as GLMTS1S2, FieldsStructured as GLMTFields



class BaseSource(object):

    def __init__(self,
                 Wavelength,
                 Polarization,
                 NA = 0.2):

        self.Wavelength = Wavelength
        self.k = 2 * np.pi / Wavelength
        self.Polarization = _Polarization(Polarization)


class MeshProperty(object):
    """Short summary.

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

    def Coupling(self, Scatterer):
        return Coupling(Scatterer = Scatterer, Detector = self)


    def Footprint(self, Scatterer, Num = 200):
        return Footprint(Scatterer = Scatterer, Detector = self, Num=Num)


    def SphericalMesh(self,
                      Sampling,
                      MaxAngle,
                      PhiOffset   = 0,
                      GammaOffset = 0,
                      Structured  = True):

        if Structured:
            x, y = np.mgrid[-50: 50: complex(Sampling), -50: 50: complex(Sampling)]
            z = 50 / np.tan(MaxAngle)
            _, phi, theta = cs.cart2sp(x, y, x*0+z)

            return phi, theta

        else:
            return FibonacciMesh(MaxAngle    = MaxAngle,
                                 Sampling    = Sampling,
                                 PhiOffset   = PhiOffset,
                                 GammaOffset = GammaOffset)

    def Plot(self):
        Name = 'Mode Field'
        ThetaMean = np.mean(self.Mesh.Theta.Degree).round(1)
        PhiMean = np.mean(self.Mesh.Phi.Degree).round(1)

        fig, (ax0, ax1) = plt.subplots(1,
                                 2,
                                 figsize=(8,4),
                                 subplot_kw = {'projection':ccrs.LambertAzimuthalEqualArea(central_latitude=PhiMean, central_longitude=ThetaMean)})

        for iter, ax in enumerate([ax0, ax1]):
            if iter == 0 : data = self.Scalar.real; Part = 'Real'
            if iter == 1 : data = self.Scalar.imag; Part = 'Imaginary'

            im = ax.tricontourf(self.Mesh.Theta.Degree,
                                    self.Mesh.Phi.Degree,
                                    data,
                                    levels=13,
                                    cmap="inferno",
                                    transform = ccrs.PlateCarree())

            ax.plot(self.Mesh.Theta.Degree,
                    self.Mesh.Phi.Degree,
                    'ko',
                    ms=0.1,
                    transform = ccrs.PlateCarree())


            plt.colorbar(mappable=im, fraction=0.046, orientation='horizontal', ax=ax)
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False)

            gl.top_labels = False
            gl.left_labels = False
            gl.right_labels = False
            gl.bottom_labels = True

            ax.set_title(f'{Part} Part {Name}')
            ax.set_ylabel(r'Angle $\phi$ [Degree]')
            ax.set_xlabel(r'Angle $\theta$ [Degree]')

        plt.show()



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

        self._Qsca, self._Qext, self._Qabs = GetEfficiencies(Index         = self.Index,
                                                             SizeParameter = self.SizeParam)


    def S1S2(self, Num=200):
        """Method compute :math:`S_1(\\phi)` & :math:`S_2(\\phi)`.
        For spherical Scatterer such as here S1 and S2 are computed as follow:

        :math:`S_1=\sum\limits_{n=1}^{n_{max}} \\frac{2n+1}{n(n+1)}(a_n\pi_n+b_n\\tau_n)`

        :math:`S_2=\sum\limits_{n=1}^{n_{max}}\\frac{2n+1}{n(n+1)}(a_n\\tau_n+b_n\pi_n)`

        Parameters
        ----------
        Num : :class:`int`
            Number of point (:math:`\\phi`) to evaluate :math:`S_1` & :math:`S_2`.

        Returns
        -------
        :class:`S1S2`
            Description of returned object.

        """

        return S1S2(Parent=self, Num=Num)



    def FarField(self, Num: int = 200):
        """Method Compute scattering Far Field.

        Fields = :math:`E_{||}(\\phi,\\theta)^2, E_{\\perp}(\\phi,\\theta)^2`.

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
        """Method Compute scattering Far Field.

        Fields = :math:`E_{||}(\\phi,\\theta)^2, E_{\\perp}(\\phi,\\theta)^2`

        The Fields are up to a constant phase value.

        :math:`\\exp{\big(-i k r \big)}`


        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\phi , \\theta`) evaluate the Fields [Num, Num].

        """

        if not self.Source.GLMT:
            return LMTFields(Index        = self.Index,
                             Diameter     = self.Diameter,
                             Wavelength   = self.Source.Wavelength,
                             nMedium      = self.nMedium,
                             Phi          = Phi,
                             Theta        = Theta,
                             Polarization = self.Source.Polarization.Radian,
                             E0           = float(self.Source.E0),
                             R            = 1.,
                             Lenght       = Phi.flatten().size)

        else:
            return GLMTFields(Index       = self.Index,
                             Diameter     = self.Diameter,
                             Wavelength   = self.Source.Wavelength,
                             nMedium      = self.nMedium,
                             Phi          = Phi,
                             Theta        = Theta,
                             Polarization = self.Source.Polarization.Radian,
                             E0           = float(self.Source.E0),
                             R            = 1.,
                             BSC          = self.Source._BSC_)




    def SPF(self, Num=100):
        """Scattering phase function.

        SPF = :math:`E_{\\parallel}(\\phi,\\theta)^2 + E_{\\perp}(\\phi,\\theta)^2`

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\theta , \\phi`) evaluate the SPF [Num, Num].

        Returns
        -------
        :class:`SPF`
            Dictionnay subclass with all pertient information as keys.

        """

        return SPF(Parent=self, Num=Num)


    def Footprint(self, Detector):
        """Method return the scattering footprint of the scatterer defined as.

        :math:`\\big| \\mathscr{F}^{-1} \\big\\{ \\tilde{ \\psi } (\\xi, \\nu),\
         \\tilde{ \\phi}_{l,m}(\\xi, \\nu)  \\big\\}(\\delta_x, \\delta_y) \\big|^2`.

        Parameters
        ----------
        Detector : type
            Description of parameter `Detector`.

        Returns
        -------
        type
            Description of returned object.

        """
        return GetFootprint(Scatterer = self, Detector = Detector)






# -
