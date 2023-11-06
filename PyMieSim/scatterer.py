#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from dataclasses import dataclass

from PyOptik import ExpData
from PyMieSim.mesh import FibonacciMesh
from PyMieSim.source import PlaneWave
from PyMieSim.tools.constants import c, epsilon0
from PyMieSim.representations import S1S2, FarField, Stokes, SPF


class GenericScatterer():
    """
    .. note::

        Generic class for scatterer

    """

    @property
    def size_parameter(self):
        return self.Bind.size_parameter

    @property
    def area(self):
        return self.Bind.area

    @property
    def Qsca(self):
        """ Scattering efficiency. """
        return self.Bind.Qsca

    @property
    def Qext(self):
        """ Extinction efficiency. """
        return self.Bind.Qext

    @property
    def Qabs(self):
        """ Absorption efficiency. """
        return self.Bind.Qabs

    @property
    def Qback(self):
        """ Backscattering efficiency. """
        return self.Bind.Qback

    @property
    def Qratio(self):
        """ Efficiency: Ratio of backscattering over total scattering. """
        return self.Bind.Qratio

    @property
    def g(self):
        """ Anisotropy factor. """
        return self.Bind.g

    @property
    def Qpr(self):
        """ Radiation pressure efficiency. """
        return self.Bind.Qpr

    @property
    def Csca(self):
        """ Scattering cross-section. """
        return (self.Bind.Csca)

    @property
    def Cext(self):
        """ Extinction cross-section. """
        return (self.Bind.Cext)

    @property
    def Cabs(self):
        """ Absorption cross-section. """
        return (self.Bind.Cabs)

    @property
    def Cpr(self):
        """ Radiation pressure cross-section. """
        return (self.Bind.Cpr)

    @property
    def Cback(self):
        """ Backscattering cross-section. """
        return (self.Bind.Cback)

    @property
    def Cratio(self):
        """ Ratio of backscattering cross-section over total scattering. """
        return (self.Bind.Cratio)

    def _FarField(self,
            phi: numpy.ndarray,
            theta: numpy.ndarray,
            r: numpy.ndarray,
             structured: bool = False) -> numpy.array:
        r"""
        .. note::

            Method Compute scattering Far Field for unstructured coordinate.

            .. math::
                \text{Fields} = E_{||}(\phi,\theta),
                                 E_{\perp}(\phi,\theta)


            The Fields are up to a constant phase value.

            .. math::
                \exp{\big(-i k r \big)}
        """

        if structured:
            return self.Bind.get_full_fields(phi.size, r=r)
        else:
            return self.Bind.get_fields(phi=phi, theta=theta, r=r)

    def get_s1s2(self, **kwargs) -> S1S2:
        r"""
        .. note::

            Method compute :math:`S_1(\phi)` and :math:`S_2(\phi)`.
            For spherical Scatterer such as here S1 and S2 are computed as follow:

            .. math::
                S_1=\sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)}(a_n \pi_n+b_n \tau_n)

                .

                S_2=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n \tau_n+b_n \pi_n)
        """

        return S1S2(parent_scatterer=self, **kwargs)

    def get_stokes(self, **kwargs) -> Stokes:
        r"""
        .. note::

            Method compute and return the Stokes parameters: I, Q, U, V.
            Those parameters are defined as:

            .. math:
                I &= \big| E_x \big|^2 + \big| E_y \big|^2

                Q &= \big| E_x \big|^2 - \big| E_y \big|^2

                U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\}

                V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\}
        """

        return Stokes(parent_scatterer=self, **kwargs)

    def get_far_field(self, **kwargs) -> FarField:
        r"""
        .. note::

            Method Compute scattering Far Field.

            .. math::
                \text{Fields} = E_{||}(\phi,\theta)^2,
                                 E_{\perp}(\phi,\theta)^2


            The Fields are up to a constant phase value:

            .. math::
                \exp{\big(-i k r \big)}
        """

        return FarField(parent_scatterer=self, **kwargs)

    def get_spf(self, **kwargs) -> SPF:
        r"""
        .. note::

            Scattering phase function.

            .. math::
                \text{SPF} = \sqrt{ E_{\parallel}(\phi,\theta)^2
                + E_{\perp}(\phi,\theta)^2 }
        """

        return SPF(parent_scatterer=self, **kwargs)

    def get_poynting_vector(self, Mesh: FibonacciMesh) -> float:
        r"""
        .. note::

            Method return the Poynting vector norm defined as:

            .. math::
                \vec{S} = \epsilon c^2 \vec{E} \times \vec{B}

        Parameters :
            Mesh : Number of voxel in the 4 pi space to compute energy flow.

        """
        Ephi, Etheta = self._FarField(
            phi=Mesh.get_phi(unit='radian'),
            theta=Mesh.get_theta(unit='radian'),
            r=1.,
            structured=False
        )

        E_norm = numpy.sqrt(numpy.abs(Ephi)**2 + numpy.abs(Etheta)**2)

        B_norm = E_norm / c

        poynting = epsilon0 * c**2 * E_norm * B_norm

        return poynting

    def get_energy_flow(self, Mesh: FibonacciMesh) -> float:
        r"""
        .. note::

            Method return energy flow defined as:

            .. math::

                W_a &= \sigma_{sca} * I_{inc}

                .

                P &= \int_{A} I dA

                .

                I &= \frac{c n \epsilon_0}{2} |E|^2

            | With:
            |     I : Energy density
            |     n  : Refractive index of the medium
            |     :math:`\epsilon_0` : Vaccum permitivity
            |     E  : Electric field
            |     \sigma_{sca}: Scattering cross section.

            More info on wikipedia link (see ref[6]).
        """

        Poynting = self.get_poynting_vector(Mesh)

        if Mesh.structured:
            Wtotal = 0.5 * numpy.sum(Poynting * Mesh.SinMesh) * Mesh.d_omega_radian

        else:
            Wtotal = 0.5 * numpy.sum(Poynting) * Mesh.d_omega_radian

        return Wtotal

    def get_cross_section(self):
        return (self.Qsca * self.area)  # similar to self.EnergyFlow(Mesh) / self.source.I

    def _assign_index_or_material(self, index, material):
        assert bool(index) ^ bool(material), logging.error("Exactly one of the parameter [index or Material] have to be assigned.")
        index = index if index is not None else material.GetRI(self.source.wavelength)
        material = material if material is not None else None
        return index, material


@dataclass()
class Sphere(GenericScatterer):
    r"""
    .. note::

        Class representing a homogeneous spherical scatterer.
    """
    diameter: float
    """ diameter of the single scatterer in unit of meter. """
    source: PlaneWave
    """ Light source object containing info on polarization and wavelength. """
    index: complex = None
    """ Refractive index of scatterer. """
    n_medium: float = 1.0
    """ Refractive index of scatterer medium. """
    material: ExpData = None
    """ Material of which the scatterer is made of. Only if index is not specified. """

    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(self.index, self.material)

        self._get_binding_()

    def _get_binding_(self) -> None:
        r"""
        .. note::

            Method call and bind c++ scatterer class

        """
        from PyMieSim.binary.SphereInterface import SPHERE

        self.Bind = SPHERE(
            wavelength=self.source.wavelength,
            amplitude=self.source.amplitude,
            diameter=self.diameter,
            index=self.index,
            n_medium=self.n_medium,
            jones_vector=self.source.linear_polarization.jones_vector.squeeze()
        )

    def an(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::

            Compute :math:`a_n` coefficient as defined in Eq:III.88 of B&B:

            .. math::
                a_n = \frac{
                \mu_{sp} \Psi_n(\alpha) \Psi_n^\prime(\beta) -
                \mu M \Psi_n^\prime(\alpha) \Psi_n(\beta)}
                {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)-
                \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)}

            With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.an()

    def bn(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::

            Compute :math:`b_n` coefficient as defined in Eq:III.89 of B&B:

            .. math::
                b_n = \frac{
                \mu M \Psi_n(\alpha) \Psi_n^\prime(\beta) -
                \mu_{sp} \Psi_n^\prime(\alpha) \Psi_n(\beta)}
                {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-
                \mu_{sp} \xi_n^\prime (\alpha) \Psi_n(\beta)}

            With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.bn()

    def cn(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::
            Compute :math:`c_n` coefficient as defined in Eq:III.90 of B&B:

            .. math::
                c_n = \frac{
                \mu_{sp} M \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) -
                \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
                {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)-
                \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)}

            With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.cn()

    def dn(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::

            Compute :math:`d_n` coefficient as defined in Eq:III.91 of B&B:

            .. math::
                d_n = \frac{
                \mu M^2 \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) -
                \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
                {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-
                \mu_{sp} M \xi_n^\prime (\alpha) \Psi_n(\beta)}

            With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        """
        return self.Bind.dn()


@dataclass()
class CoreShell(GenericScatterer):
    r"""
    .. note::
        Class representing a core+shell spherical scatterer.

    """

    core_diameter: float
    """ diameter of the core of the single scatterer [m]. """
    shell_width: float
    """ diameter of the shell of the single scatterer [m]. """
    source: PlaneWave
    """ Light source object containing info on polarization and wavelength. """
    core_index: complex = None
    """ Refractive index of the core of the scatterer. """
    shell_index: complex = None
    """ Refractive index of the shell of the scatterer. """
    core_material: ExpData = None
    """ Core material of which the scatterer is made of. Only if core_index is not specified.  """
    shell_material: ExpData = None
    """ Shell material of which the scatterer is made of. Only if shell_index is not specified.  """
    n_medium: float = 1.0
    """ Refractive index of scatterer medium. """

    def __post_init__(self):
        self.core_index, self.core_material = self._assign_index_or_material(self.core_index, self.core_material)

        self.shell_index, self.shell_material = self._assign_index_or_material(self.shell_index, self.shell_material)

        self.shell_diameter = self.core_diameter + self.shell_width

        self._get_binding_()

    def _get_binding_(self) -> None:
        r"""
        .. note::
            Method call and bind c++ scatterer class
        """
        from PyMieSim.binary.CoreShellInterface import CORESHELL

        self.Bind = CORESHELL(
            shell_index=self.shell_index,
            core_index=self.core_index,
            shell_width=self.shell_width,
            core_diameter=self.core_diameter,
            wavelength=self.source.wavelength,
            n_medium=self.n_medium,
            jones_vector=self.source.linear_polarization.jones_vector.squeeze(),
            amplitude=self.source.amplitude
        )

    def an(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::
            Compute :math:`a_n` coefficient

        """
        if max_order is None:
            return self.Bind.an()
        else:
            return self.Bind._an(max_order)

    def bn(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::
            Compute :math:`b_n` coefficient.

        """
        if max_order is None:
            return self.Bind.bn()
        else:
            return self.Bind._bn(max_order)


@dataclass()
class Cylinder(GenericScatterer):
    r"""
    .. note::
        Class representing a right angle cylindrical scatterer.
    """

    diameter: float
    """ diameter of the single scatterer in unit of meter. """
    source: PlaneWave
    """ Light source object containing info on polarization and wavelength. """
    index: complex = None
    """ Refractive index of scatterer. """
    n_medium: float = 1.0
    """ Material of which the scatterer is made of. Only if index is not specified. """
    material: ExpData = None
    """ Refractive index of scatterer medium. """

    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(self.index, self.material)

        self._get_binding_()

    def _get_binding_(self) -> None:
        r"""
        .. note::

            Method call and bind c++ scatterer class

        """
        from PyMieSim.binary.CylinderInterface import CYLINDER

        self.Bind = CYLINDER(
            index=self.index,
            diameter=self.diameter,
            wavelength=self.source.wavelength,
            n_medium=self.n_medium,
            amplitude=self.source.amplitude,
            jones_vector=self.source.linear_polarization.jones_vector.squeeze()
        )

    def a1n(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::
            Compute :math:`a_n` coefficient as defined ref[5]:

            .. math::
                a_n = \frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
                { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

            | With :math:`m` being the refractive index of the medium and
            |      :math:`m_t` being the refractive index of the index.

        """

        return self.Bind.a1n()

    def a2n(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::
            Compute :math:`a_n` coefficient as defined ref[5]:

            .. math::
                a_n = \frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
                { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

            | With :math:`m` being the refractive index of the medium and
            |      :math:`m_t` being the refractive index of the index.

        """

        return self.Bind.a2n()

    def b1n(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::
            Compute :math:`b_n` coefficient as defined in ref[5]:

            .. math::
                b_n = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
                { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

            | With :math:`m` being the refractive index of the medium and
            |      :math:`m_t` being the refractive index of the index.

        """

        return self.Bind.b1n()

    def b2n(self, max_order: int = None) -> numpy.array:
        r"""
        .. note::
            Compute :math:`b_n` coefficient as defined in ref[5]:

            .. math::
                b_n = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
                { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

            | With :math:`m` being the refractive index of the medium and
            |      :math:`m_t` being the refractive index of the index.

        """

        return self.Bind.b2n()

# -
