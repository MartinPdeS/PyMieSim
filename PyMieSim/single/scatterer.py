#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.base_class import BaseMaterial

import numpy
from tabulate import tabulate
from typing import Union, Optional, Any  # Any is for complex as no

from PyMieSim.single.representations import S1S2, FarField, Stokes, SPF, Footprint
from PyMieSim.single import source  # noqa: F401
from pydantic.dataclasses import dataclass


config_dict = dict(
    arbitrary_types_allowed=True,
    kw_only=True,
    slots=True,
    extra='forbid'
)


class GenericScatterer:
    """
    A generic class for a scatterer, providing properties and methods to compute various scattering-related quantities.
    """

    def print_properties(self) -> None:
        """Prints a table of the scatterer's properties and their values."""
        data = [getattr(self, name) for name in self.property_names]
        property_dict = {"Property": self.property_names, "Value": data}

        table = tabulate(property_dict, headers="keys")
        print(table)

    @property
    def size_parameter(self) -> float:
        """Returns the size parameter of the scatterer."""
        return self.binding.size_parameter

    @property
    def area(self) -> float:
        """Returns the area of the scatterer."""
        return self.binding.area

    @property
    def Qsca(self) -> float:
        """Returns the scattering efficiency."""
        return self.binding.Qsca

    @property
    def Qext(self) -> float:
        """Returns the extinction efficiency."""
        return self.binding.Qext

    @property
    def Qabs(self) -> float:
        """Returns the absorption efficiency."""
        return self.binding.Qabs

    @property
    def Qback(self) -> float:
        """Returns the backscattering efficiency."""
        return self.binding.Qback

    @property
    def Qforward(self) -> float:
        """Returns the forward-scattering efficiency."""
        return self.binding.Qforward

    @property
    def Qratio(self) -> float:
        """Returns the efficiency ratio of backscattering over total scattering."""
        return self.binding.Qratio

    @property
    def g(self) -> float:
        """Returns the anisotropy factor."""
        return self.binding.g

    @property
    def Qpr(self) -> float:
        """Returns the radiation pressure efficiency."""
        return self.binding.Qpr

    @property
    def Csca(self) -> float:
        """Returns the scattering cross-section."""
        return self.binding.Csca

    @property
    def Cext(self) -> float:
        """Returns the extinction cross-section."""
        return self.binding.Cext

    @property
    def Cabs(self) -> float:
        """Returns the absorption cross-section."""
        return self.binding.Cabs

    @property
    def Cpr(self) -> float:
        """Returns the radiation pressure cross-section."""
        return self.binding.Cpr

    @property
    def Cback(self) -> float:
        """Returns the backscattering cross-section."""
        return self.binding.Cback

    @property
    def Cforward(self) -> float:
        """Returns the forward-scattering cross-section."""
        return self.binding.Cforward

    @property
    def Cratio(self) -> float:
        """Returns the ratio of backscattering cross-section over total scattering."""
        return self.binding.Cratio

    def get_farfields_array(self, phi: numpy.ndarray, theta: numpy.ndarray, r: numpy.ndarray) -> numpy.ndarray:
        """
        Computes the scattering far field for unstructured coordinates.

        The method computes the fields up to a constant phase value.

        Args:
            phi (numpy.ndarray): The phi angles in radians.
            theta (numpy.ndarray): The theta angles in radians.
            r (numpy.ndarray): The radial distances.

        Returns:
            numpy.ndarray: The computed far fields.
        """
        return self.binding.get_fields(phi=phi, theta=theta, r=r)

    def get_s1s2(self, **kwargs) -> S1S2:
        r"""
        Computes the S1 and S2 parameters.

        These parameters are computed for spherical scatterers using the following formulas:

        .. math::
            S_1=\sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)}(a_n \pi_n+b_n \tau_n) \\
            S_2=\sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)}(a_n \tau_n+b_n \pi_n)

        Args:
            **kwargs: Additional keyword arguments.

        Returns:
            S1S2: The computed S1 and S2 parameters.
        """
        return S1S2(scatterer=self, **kwargs)

    def get_stokes(self, **kwargs) -> Stokes:
        r"""
        Computes the four Stokes components: I, Q, U, and V.

        These components are defined as:

        .. math::
            I &= \big| E_x \big|^2 + \big| E_y \big|^2 \\
            Q &= \big| E_x \big|^2 - \big| E_y \big|^2 \\
            U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\} \\
            V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\}

        Args:
            **kwargs: Additional keyword arguments.

        Returns:
            Stokes: The computed Stokes parameters.
        """
        return Stokes(scatterer=self, **kwargs)

    def get_far_field(self, **kwargs) -> FarField:
        r"""
        Computes the scattering far fields.

        The far fields are computed as:

        .. math::
            \text{Fields} = E_{||}(\phi,\theta)^2, E_{\perp}(\phi,\theta)^2

        The fields are up to a constant phase value:

        .. math::
            \exp{\big(-i k r \big)}

        Args:
            **kwargs: Additional keyword arguments.

        Returns:
            FarField: The computed far fields.
        """
        return FarField(scatterer=self, **kwargs)

    def get_spf(self, **kwargs) -> SPF:
        r"""
        Computes the scattering phase function.

        The scattering phase function is computed as:

        .. math::
            \text{SPF} = \sqrt{ E_{\parallel}(\phi,\theta)^2 + E_{\perp}(\phi,\theta)^2 }

        Args:
            **kwargs: Additional keyword arguments.

        Returns:
            SPF: The computed scattering phase function.
        """
        return SPF(scatterer=self, **kwargs)

    def get_footprint(self, detector) -> Footprint:
        r"""
        Computes the footprint of the scattered light coupling with the detector.

        The footprint is computed as:

        .. math::
            \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi } (\xi, \nu),\
                   \tilde{ \phi}_{l,m}(\xi, \nu)  \big\}
                   (\delta_x, \delta_y) \big|^2

        Args:
            detector (GenericDetector): The detector object.

        Returns:
            Footprint: The computed scatterer footprint.
        """
        return Footprint(scatterer=self, detector=detector)

    def get_cross_section(self) -> float:
        """
        Computes the scattering cross-section.

        Returns:
            float: The scattering cross-section, calculated as `Qsca * area`.
        """
        return self.Qsca * self.area

    def _assign_index_or_material(self, index: Optional[Any], material: Optional[BaseMaterial]) -> tuple:
        """
        Assigns the refractive index or material.

        Either `index` or `material` must be specified, but not both.

        Args:
            index (Optional[Any]): The refractive index.
            material (Optional[BaseMaterial]): The material.

        Returns:
            tuple: A tuple containing the refractive index and material.

        Raises:
            ValueError: If neither or both of `index` and `material` are specified.
        """
        if (index is None) == (material is None):
            raise ValueError("Either index or material must be specified, but not both.")

        if index is None:
            index = material.compute_refractive_index(self.source.wavelength)

        if not numpy.isscalar(index) and len(index) == 1:
            index = index[0]

        return index, material


@dataclass(config=config_dict)
class Sphere(GenericScatterer):
    """
    Class representing a homogeneous spherical scatterer.

    Attributes:
        diameter (float): Diameter of the single scatterer in unit of meter.
        source (Union[source.PlaneWave, source.Gaussian]): Light source object containing info on polarization and wavelength.
        index (Optional[Any]): Refractive index of scatterer. Default is None.
        medium_index (float): Refractive index of scatterer medium. Default is 1.0.
        material (Optional[BaseMaterial]): BaseMaterial of which the scatterer is made, if index is not specified. Default is None.
    """
    diameter: float
    source: Union[source.PlaneWave, source.Gaussian]
    index: Optional[Any] = None
    medium_index: Optional[float] = None
    medium_material: Optional[BaseMaterial] = None
    material: Optional[BaseMaterial] = None

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qpr",
        "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cpr"
    ]
    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(self.index, self.material)

        self.medium_index, self.medium_material = self._assign_index_or_material(self.medium_index, self.medium_material)
        self.set_binding()

    def set_binding(self) -> None:
        """
        Binds the Python representation of the sphere to its C++ counterpart using provided properties.
        """
        from PyMieSim.binary.SphereInterface import SPHERE

        self.binding = SPHERE(
            diameter=self.diameter,
            index=self.index,
            medium_index=self.medium_index,
            source=self.source.binding
        )

    def an(self, max_order: Optional[int] = 0) -> numpy.ndarray:
        r"""
        Compute :math:`a_n` coefficient as defined in Eq:III.88 of B&B:

        .. math::
            a_n = \frac{
            \mu_{sp} \Psi_n(\alpha) \Psi_n^\prime(\beta) -
            \mu M \Psi_n^\prime(\alpha) \Psi_n(\beta)}
            {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)-
            \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)}

        With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`a_n` coefficients.
        """
        return self.binding.an(max_order)

    def bn(self, max_order: Optional[int] = 0) -> numpy.ndarray:
        r"""
        Compute :math:`b_n` coefficient as defined in Eq:III.89 of B&B:

        .. math::
            b_n = \frac{
            \mu M \Psi_n(\alpha) \Psi_n^\prime(\beta) -
            \mu_{sp} \Psi_n^\prime(\alpha) \Psi_n(\beta)}
            {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-
            \mu_{sp} \xi_n^\prime (\alpha) \Psi_n(\beta)}

        With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`b_n` coefficients.
        """
        return self.binding.bn(max_order)

    def cn(self, max_order: Optional[int] = 0) -> numpy.ndarray:
        r"""
        For future purpose only!
        Compute :math:`c_n` coefficient as defined in Eq:III.90 of B&B:

        .. math::
            c_n = \frac{
            \mu_{sp} M \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) -
            \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
            {\mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta)-
            \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)}

        With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`c_n` coefficients.
        """
        return self.binding.cn(max_order)

    def dn(self, max_order: Optional[int] = 0) -> numpy.ndarray:
        r"""
        For future purpose only!
        Compute :math:`d_n` coefficient as defined in Eq:III.91 of B&B:

        .. math::
            d_n = \frac{
            \mu M^2 \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) -
            \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]}
            {\mu M \xi_n(\alpha) \Psi_n^\prime(\beta)-
            \mu_{sp} M \xi_n^\prime (\alpha) \Psi_n(\beta)}

        With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`d_n` coefficients.
        """
        return self.binding.dn(max_order)


@dataclass(config=config_dict)
class CoreShell(GenericScatterer):
    """
    Class representing a core/shell spherical scatterer.

    Attributes:
        core_diameter (float): Diameter of the core of the single scatterer [m].
        shell_width (float): Diameter of the shell of the single scatterer [m].
        source (Union[source.PlaneWave, source.Gaussian]): Light source object containing info on polarization and wavelength.
        core_index (Optional[Any]): Refractive index of the core of the scatterer. Default is None.
        shell_index (Optional[Any]): Refractive index of the shell of the scatterer. Default is None.
        core_material (Optional[BaseMaterial]): Core material of which the scatterer is made of, if core_index is not specified. Default is None.
        shell_material (Optional[BaseMaterial]): Shell material of which the scatterer is made of, if shell_index is not specified. Default is None.
        medium_index (float): Refractive index of the scatterer medium. Default is 1.0.
    """

    core_diameter: float
    shell_width: float
    source: Union[source.PlaneWave, source.Gaussian]
    core_index: Optional[Any] = None
    shell_index: Optional[Any] = None
    core_material: Optional[BaseMaterial] = None
    shell_material: Optional[BaseMaterial] = None
    medium_index: Optional[float] = None
    medium_material: Optional[BaseMaterial] = None

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qpr",
        "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cpr"
    ]

    def __post_init__(self):
        self.core_index, self.core_material = self._assign_index_or_material(self.core_index, self.core_material)
        self.shell_index, self.shell_material = self._assign_index_or_material(self.shell_index, self.shell_material)
        self.medium_index, self.medium_material = self._assign_index_or_material(self.medium_index, self.medium_material)
        self.shell_diameter = self.core_diameter + self.shell_width
        self.set_binding()

    def set_binding(self) -> None:
        """
        Bind the C++ scatterer class.
        """
        from PyMieSim.binary.CoreShellInterface import CORESHELL

        self.binding = CORESHELL(
            shell_index=self.shell_index,
            core_index=self.core_index,
            shell_width=self.shell_width,
            core_diameter=self.core_diameter,
            medium_index=self.medium_index,
            source=self.source.binding
        )

    def an(self, max_order: Optional[int] = 0) -> numpy.ndarray:
        r"""
        Compute :math:`a_n` coefficient.

        If max_order is set to zero, the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`a_n` coefficients.
        """
        return self.binding.an(max_order)

    def bn(self, max_order: Optional[int] = 0) -> numpy.ndarray:
        r"""
        Compute :math:`b_n` coefficient.

        If max_order is set to zero, the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`b_n` coefficients.
        """
        return self.binding.bn(max_order)


@dataclass(config=config_dict)
class Cylinder(GenericScatterer):
    """
    Class representing a right angle cylindrical scatterer.

    Attributes:
        diameter (float): Diameter of the single scatterer in unit of meter.
        source (Union[source.PlaneWave, source.Gaussian]): Light source object containing info on polarization and wavelength.
        index (Optional[Any]): Refractive index of scatterer. Default is None.
        medium_index (float): Refractive index of scatterer medium. Default is 1.0.
        material (Union[BaseMaterial, None]): BaseMaterial of which the scatterer is made, if index is not specified. Default is None.
    """

    diameter: float
    source: Union[source.PlaneWave, source.Gaussian]
    index: Optional[Any] = None
    medium_index: Optional[float] = None
    medium_material: Optional[BaseMaterial] = None
    material: Union[BaseMaterial, None] = None

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs",
        "Csca", "Cext", "Cabs"
    ]

    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(index=self.index, material=self.material)
        self.medium_index, self.medium_material = self._assign_index_or_material(index=self.medium_index, material=self.medium_material)
        self.set_binding()

    def set_binding(self) -> None:
        """
        Binds the Python representation of the cylinder to its C++ counterpart using provided properties.
        """
        from PyMieSim.binary.CylinderInterface import CYLINDER

        self.binding = CYLINDER(
            index=self.index,
            diameter=self.diameter,
            medium_index=self.medium_index,
            source=self.source.binding
        )

    def a1n(self, max_order: Optional[int] = 0) -> numpy.array:
        r"""
        Compute :math:`a_{1n}` coefficient as defined in ref[5]:

        .. math::
            a_{1n} = \frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
                  { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

        With :math:`m` being the refractive index of the medium and
        :math:`m_t` being the refractive index of the scatterer.

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`a_{1n}` coefficients.
        """
        return self.binding.a1n(max_order)

    def a2n(self, max_order: Optional[int] = 0) -> numpy.array:
        r"""
        Compute :math:`a_{2n}` coefficient as defined in ref[5]:

        .. math::
            a_{2n} = \frac{ m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x) }
                  { m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x) }

        With :math:`m` being the refractive index of the medium and
        :math:`m_t` being the refractive index of the scatterer.

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`a_{2n}` coefficients.
        """
        return self.binding.a2n(max_order)

    def b1n(self, max_order: Optional[int] = 0) -> numpy.array:
        r"""
        Compute :math:`b_{1n}` coefficient as defined in ref[5]:

        .. math::
            b_{1n} = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
                  { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

        With :math:`m` being the refractive index of the medium and
        :math:`m_t` being the refractive index of the scatterer.

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`b_{1n}` coefficients.
        """
        return self.binding.b1n(max_order)

    def b2n(self, max_order: Optional[int] = 0) -> numpy.array:
        r"""
        Compute :math:`b_{2n}` coefficient as defined in ref[5]:

        .. math::
            b_{2n} = \frac{ m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x) }
                  { m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x) }

        With :math:`m` being the refractive index of the medium and
        :math:`m_t` being the refractive index of the scatterer.

        If max_order is set to zero then the maximum order output is calculated using the Wiscombe criterion.

        Args:
            max_order (Optional[int]): The maximum order of the coefficient. Default is 0.

        Returns:
            numpy.ndarray: Array of :math:`b_{2n}` coefficients.
        """
        return self.binding.b2n(max_order)

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

# -
