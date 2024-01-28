#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyOptik import DataMeasurement, Sellmeier

import numpy
import logging
from dataclasses import dataclass
from tabulate import tabulate


from PyMieSim.mesh import FibonacciMesh
from PyMieSim.single.source import PlaneWave, Gaussian
from PyMieSim.single.representations import S1S2, FarField, Stokes, SPF, Footprint

c = 299792458.0  #: Speed of light in vacuum (m/s).
epsilon0 = 8.854187817620389e-12  #: Vacuum permittivity (F/m).


class GenericScatterer():
    """
    Generic class for scatterer
    """
    def print_properties(self) -> None:
        property_names = [
            "size_parameter",
            "area",
            "index",
            "Qsca",
            "Qext",
            "Qabs",
            "Qback",
            "Qratio",
            "Qpr",
            "Csca",
            "Cext",
            "Cabs",
            "Cback",
            "Cratio",
            "Cpr",
            "g",
        ]

        data = [getattr(self, name) for name in property_names]
        property_dict = {"Property": property_names, "value": data}

        table = tabulate(
            property_dict,
            headers="keys"
        )
        print(table)

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
        return self.Bind.Qback / self.Bind.Qsca

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
        return self.Bind.Csca

    @property
    def Cext(self):
        """ Extinction cross-section. """
        return self.Bind.Cext

    @property
    def Cabs(self):
        """ Absorption cross-section. """
        return self.Bind.Cabs

    @property
    def Cpr(self):
        """ Radiation pressure cross-section. """
        return self.Bind.Cpr

    @property
    def Cback(self):
        """ Backscattering cross-section. """
        return self.Bind.Cback

    @property
    def Cratio(self):
        """ Ratio of backscattering cross-section over total scattering. """
        return self.Bind.Cback / self.Bind.Csca

    def get_farfields_array(
            self,
            phi: numpy.ndarray,
            theta: numpy.ndarray,
            r: numpy.ndarray,
            structured: bool = False) -> numpy.array:
        r"""
        Method Compute scattering Far Field for unstructured coordinate.

        .. math::
            \text{Fields} = E_{||}(\phi,\theta), E_{\perp}(\phi,\theta)


        The Fields are up to a constant phase value.

        .. math::
            \exp{\big(-i k r \big)}

        :param      phi:         The phi array
        :type       phi:         numpy.ndarray
        :param      theta:       The theta array
        :type       theta:       numpy.ndarray
        :param      r:           The radial array
        :type       r:           numpy.ndarray
        :param      structured:  Indicates if computing mesh is structured or not
        :type       structured:  bool

        :returns:   The far fields
        :rtype:     numpy.ndarray
        """
        if structured:
            return self.Bind.get_full_fields(phi.size, r=r)
        else:
            return self.Bind.get_fields(phi=phi, theta=theta, r=r)

    def get_s1s2(self, **kwargs) -> S1S2:
        r"""
        Method compute :math:`S_1(\phi)` and :math:`S_2(\phi)`.
        For spherical Scatterer such as here S1 and S2 are computed as follow:

        .. math::
            S_1=\sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)}(a_n \pi_n+b_n \tau_n) \\[10pt]

            S_2=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n \tau_n+b_n \pi_n) \\[10pt]

        :param      kwargs:  The keywords arguments
        :type       kwargs:  dictionary

        :returns:   The S1 and S2 parameters
        :rtype:     S1S2
        """
        return S1S2(scatterer=self, **kwargs)

    def get_stokes(self, **kwargs) -> Stokes:
        r"""
        Returns the four Stokes components. The method compute the Stokes parameters: I, Q, U, V.
        Those parameters are defined as:

        .. math:
            I &= \big| E_x \big|^2 + \big| E_y \big|^2 \\[10pt]

            Q &= \big| E_x \big|^2 - \big| E_y \big|^2 \\[10pt]

            U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\} \\[10pt]

            V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\} \\[10pt]

        :param      kwargs:  The keywords arguments
        :type       kwargs:  dictionary

        :returns:   The stokes.
        :rtype:     Stokes
        """
        return Stokes(scatterer=self, **kwargs)

    def get_far_field(self, **kwargs) -> FarField:
        r"""
        Returns the scattering far-fields defined as.

        .. math::
            \text{Fields} = E_{||}(\phi,\theta)^2, E_{\perp}(\phi,\theta)^2


        The Fields are up to a constant phase value:

        .. math::
            \exp{\big(-i k r \big)}

        :param      kwargs:  The keywords arguments
        :type       kwargs:  dictionary

        :returns:   The far field.
        :rtype:     FarField
        """
        return FarField(scatterer=self, **kwargs)

    def get_spf(self, **kwargs) -> SPF:
        r"""
        Returns the scattering phase function.

        .. math::
            \text{SPF} = \sqrt{ E_{\parallel}(\phi,\theta)^2
            + E_{\perp}(\phi,\theta)^2 }

        :param      kwargs:  The keywords arguments
        :type       kwargs:  dictionary

        :returns:   The scattering phase function.
        :rtype:     SPF
        """
        return SPF(scatterer=self, **kwargs)

    def get_footprint(self, detector) -> Footprint:
        r"""
        Return the footprint of the scattererd light coupling with the
        detector as computed as:

        .. math::
            \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi } (\xi, \nu),\
                   \tilde{ \phi}_{l,m}(\xi, \nu)  \big\}
                   (\delta_x, \delta_y) \big|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        :param      detector:  The detector
        :type       detector:  GenericDetector

        :returns:   The scatterer footprint.
        :rtype:     Footprint
        """
        return Footprint(scatterer=self, detector=detector)

    def get_poynting_vector(self, mesh: FibonacciMesh) -> float:
        r"""

        Method return the Poynting vector norm defined as:

        .. math::
            \vec{S} = \epsilon c^2 \vec{E} \times \vec{B}

        Parameters :
            Mesh : Number of voxel in the 4 pi space to compute energy flow.

        """
        Ephi, Etheta = self.get_farfields_array(
            phi=mesh.phi,
            theta=mesh.theta,
            r=1.,
            structured=False
        )

        E_norm = numpy.sqrt(numpy.abs(Ephi)**2 + numpy.abs(Etheta)**2)

        B_norm = E_norm / c

        poynting = epsilon0 * c**2 * E_norm * B_norm

        return poynting

    def get_energy_flow(self, mesh: FibonacciMesh) -> float:
        r"""
        Returns energy flow defined as:

        .. math::
            W_a &= \sigma_{sca} * I_{inc} \\[10pt]
            P &= \int_{A} I dA \\[10pt]
            I &= \frac{c n \epsilon_0}{2} |E|^2 \\[10pt]

        | With:
        |     I : Energy density
        |     n  : Refractive index of the medium
        |     :math:`\epsilon_0` : Vaccum permitivity
        |     E  : Electric field
        |     \sigma_{sca}: Scattering cross section.

        More info on wikipedia link (see ref[6]).

        :param      Mesh:  The mesh
        :type       Mesh:  FibonacciMesh

        :returns:   The energy flow.
        :rtype:     float
        """
        Poynting = self.get_poynting_vector(mesh)

        total_power = 0.5 * numpy.sum(Poynting) * mesh.d_omega

        return total_power

    def get_cross_section(self):
        return (self.Qsca * self.area)  # similar to self.EnergyFlow(Mesh) / self.source.I

    def _assign_index_or_material(self, index, material) -> tuple:
        assert bool(index) ^ bool(material), logging.error("Exactly one of the parameter [index or Material] have to be assigned.")
        index = index if index is not None else material.get_refractive_index(self.source.wavelength)
        material = material if material is not None else None

        if not numpy.isscalar(index) and len(index) == 1:
            return index[0], material

        return index, material


@dataclass()
class Sphere(GenericScatterer):
    """ Class representing a homogeneous spherical scatterer """
    diameter: float
    """ diameter of the single scatterer in unit of meter. """
    source: PlaneWave | Gaussian
    """ Light source object containing info on polarization and wavelength. """
    index: complex = None
    """ Refractive index of scatterer. """
    n_medium: float = 1.0
    """ Refractive index of scatterer medium. """
    material: DataMeasurement | Sellmeier = None
    """ Material of which the scatterer is made of. Only if index is not specified. """

    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(self.index, self.material)

        self.set_cpp_binding()

    def set_cpp_binding(self) -> None:
        """
        Method call and bind c++ scatterer class

        """
        from PyMieSim.binary.SphereInterface import SPHERE

        self.Bind = SPHERE(
            wavelength=self.source.wavelength,
            amplitude=self.source.amplitude,
            diameter=self.diameter,
            index=self.index,
            n_medium=self.n_medium,
            jones_vector=self.source.polarization.jones_vector.squeeze(),
        )

    def an(self, max_order: int = None) -> numpy.array:
        r"""
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
    """
    Class representing a core/shell spherical scatterer.
    """

    core_diameter: float
    """ diameter of the core of the single scatterer [m]. """
    shell_width: float
    """ diameter of the shell of the single scatterer [m]. """
    source: PlaneWave | Gaussian
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

        self.set_cpp_binding()

    def set_cpp_binding(self) -> None:
        """
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
            jones_vector=self.source.polarization.jones_vector.squeeze(),
            amplitude=self.source.amplitude
        )

    def an(self, max_order: int = None) -> numpy.array:
        r"""
        Compute :math:`a_n` coefficient
        """
        if max_order is None:
            return self.Bind.an()
        else:
            return self.Bind._an(max_order)

    def bn(self, max_order: int = None) -> numpy.array:
        r"""
        Compute :math:`b_n` coefficient.
        """
        if max_order is None:
            return self.Bind.bn()
        else:
            return self.Bind._bn(max_order)


@dataclass()
class Cylinder(GenericScatterer):
    """
    Class representing a right angle cylindrical scatterer.
    """

    diameter: float
    """ diameter of the single scatterer in unit of meter. """
    source: PlaneWave | Gaussian
    """ Light source object containing info on polarization and wavelength. """
    index: complex = None
    """ Refractive index of scatterer. """
    n_medium: float = 1.0
    """ Material of which the scatterer is made of. Only if index is not specified. """
    material: ExpData = None
    """ Refractive index of scatterer medium. """

    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(
            index=self.index,
            material=self.material
        )

        self.set_cpp_binding()

    def set_cpp_binding(self) -> None:
        """
        Method call and bind c++ scatterer class

        """
        from PyMieSim.binary.CylinderInterface import CYLINDER

        self.Bind = CYLINDER(
            index=self.index,
            diameter=self.diameter,
            wavelength=self.source.wavelength,
            n_medium=self.n_medium,
            amplitude=self.source.amplitude,
            jones_vector=self.source.polarization.jones_vector.squeeze()
        )

    def a1n(self, max_order: int = None) -> numpy.array:
        r"""
        Compute :math:`a_n` coefficient as defined ref[5]:

        .. math::
            a_n = \frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
            { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

        | With :math:`m` being the refractive index of the medium and
        |      :math:`m_t` being the refractive index of the index.

        :param      max_order:  The maximum order
        :type       max_order:  int

        :returns:   The first electric mutlipole amplitude
        :rtype:     numpy.array
        """
        return self.Bind.a1n()

    def a2n(self, max_order: int = None) -> numpy.array:
        r"""
        Compute :math:`a_n` coefficient as defined ref[5]:

        .. math::
            a_n = \frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
            { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

        | With :math:`m` being the refractive index of the medium and
        |      :math:`m_t` being the refractive index of the index.

        :param      max_order:  The maximum order
        :type       max_order:  int

        :returns:   The second electric mutlipole amplitude
        :rtype:     numpy.array
        """
        return self.Bind.a2n()

    def b1n(self, max_order: int = None) -> numpy.array:
        r"""
        Compute :math:`b_n` coefficient as defined in ref[5]:

        .. math::
            b_n = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
            { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

        | With :math:`m` being the refractive index of the medium and
        |      :math:`m_t` being the refractive index of the index.


        :param      max_order:  The maximum order
        :type       max_order:  int

        :returns:   The first magnetic mutlipole amplitude
        :rtype:     numpy.array
        """
        return self.Bind.b1n()

    def b2n(self, max_order: int = None) -> numpy.array:
        r"""
        Compute :math:`b_n` coefficient as defined in ref[5]:

        .. math::
            b_n = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
            { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

        | With :math:`m` being the refractive index of the medium and
        |      :math:`m_t` being the refractive index of the index.

        :param      max_order:  The maximum order
        :type       max_order:  int

        :returns:   The second magnetic mutlipole amplitude
        :rtype:     numpy.array
        """
        return self.Bind.b2n()

    @property
    def Cback(self) -> None:
        raise NotImplementedError

    @property
    def Qback(self) -> None:
        raise NotImplementedError

    @property
    def Cratio(self) -> None:
        raise NotImplementedError

    @property
    def Qratio(self) -> None:
        raise NotImplementedError

    @property
    def Cpr(self) -> None:
        raise NotImplementedError

    @property
    def Qpr(self) -> None:
        raise NotImplementedError

# -
