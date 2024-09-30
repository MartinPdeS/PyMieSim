#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.base_class import BaseMaterial

import numpy
from typing import Optional, Any  # Any is for complex as no
from PyMieSim.single.source.base import BaseSource
from pydantic.dataclasses import dataclass
from pydantic import field_validator
from PyMieSim.units import Quantity, meter, RIU
from PyMieSim.single.scatterer.base import GenericScatterer, config_dict



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
    diameter: Quantity
    source: BaseSource
    index: Optional[Quantity] = None
    medium_index: Optional[Quantity] = None
    medium_material: Optional[BaseMaterial] = None
    material: Optional[BaseMaterial] = None

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qpr",
        "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cpr"
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
        self.index, self.material = self._assign_index_or_material(self.index, self.material)

        self.medium_index, self.medium_material = self._assign_index_or_material(self.medium_index, self.medium_material)

        self.set_binding()


    def set_binding(self) -> None:
        """
        Binds the Python representation of the sphere to its C++ counterpart using provided properties.
        """
        from PyMieSim.binary.SphereInterface import SPHERE

        self.binding = SPHERE(
            diameter=self.diameter.to_base_units().magnitude,
            index=self.index.to_base_units().magnitude,
            medium_index=self.medium_index.to_base_units().magnitude,
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