#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.base_class import BaseMaterial

import numpy
from tabulate import tabulate
from typing import Optional, Any  # Any is for complex as no

from PyMieSim.single.representations import S1S2, FarField, Stokes, SPF, Footprint
from PyMieSim.units import RIU

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
            kwargs Additional keyword arguments.

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
            kwargs Additional keyword arguments.

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
            kwargs Additional keyword arguments.

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
            kwargs Additional keyword arguments.

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
            index = material.compute_refractive_index(
                self.source.wavelength.to_base_units().magnitude
            )

            from typing import Iterable
            if isinstance(index, Iterable):
                index = index[0]

            index = index * RIU

        return index, material