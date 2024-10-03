#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.base_class import BaseMaterial

import numpy
from typing import Optional
from pydantic.dataclasses import dataclass
from pydantic import field_validator
from PyMieSim.units import Quantity, meter, RIU
from PyMieSim.single.scatterer.base import BaseScatterer, config_dict



@dataclass(config=config_dict, kw_only=True)
class Cylinder(BaseScatterer):
    """
    Represents a cylindrical scatterer used for scattering simulations in optical systems.

    Parameters
    ----------
        diameter (Quantity): Diameter of the cylindrical scatterer, given in meters.
        property (Quantity | BaseMaterial): Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer. Only one of these should be provided at a time to specify the core characteristic.

    Attributes
    ----------
        size_parameter: Dimensionless size parameter of the scatterer.
        area: Effective scattering area.
        g: Asymmetry factor, representing the average cosine of the scattering angle.
        Qsca: Scattering efficiency, the ratio of the scattering cross-section to the geometric cross-section.
        Qext: Extinction efficiency, the sum of scattering and absorption efficiencies.
        Qabs: Absorption efficiency, quantifying the portion of light absorbed by the scatterer.
        Csca: Scattering cross-section, the effective area that scatters light.
        Cext: Extinction cross-section, the effective area accounting for both scattering and absorption.
        Cabs: Absorption cross-section, representing the area that absorbs light.
    """


    diameter: Quantity
    property: Quantity | BaseMaterial

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs",
        "Csca", "Cext", "Cabs"
    ]

    @field_validator('diameter', mode='before')
    def validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units.
        """
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return value

    @field_validator('index', 'medium_index', mode='before')
    def validate_riu_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with RIU units.
        """
        if value is None:
            return None

        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(RIU):
            raise ValueError(f"{value} must have RIU units.")

        return value

    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(self.property)

        super().__post_init__()

    def set_binding(self) -> None:
        """
        Binds the Python representation of the cylinder to its C++ counterpart using provided properties.
        """
        from PyMieSim.binary.CylinderInterface import CYLINDER

        self.binding = CYLINDER(
            index=self.index.to_base_units().magnitude,
            diameter=self.diameter.to_base_units().magnitude,
            medium_index=self.medium_index.to_base_units().magnitude,
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