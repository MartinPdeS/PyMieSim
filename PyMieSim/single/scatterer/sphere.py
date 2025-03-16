#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.material.base_class import BaseMaterial

import numpy
from typing import Optional
from pydantic.dataclasses import dataclass
from PyMieSim.units import Quantity
from PyMieSim.single.scatterer.base import BaseScatterer, config_dict


@dataclass(config=config_dict, kw_only=True)
class Sphere(BaseScatterer):
    """
    Class representing a homogeneous spherical scatterer.

    Parameters
    ----------
    diameter : Quantity
        Diameter of the cylindrical scatterer, given in meters.
    property : Quantity or BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer. Only one can be provided.

    """
    diameter: Quantity
    property: Quantity | BaseMaterial

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qpr",
        "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cpr"
    ]

    def __post_init__(self):
        self.index, self.material = self._assign_index_or_material(self.property)

        self.cross_section = numpy.pi * (self.diameter / 2) ** 2

        super().__post_init__()

    def set_binding(self) -> None:
        """
        Binds the Python representation of the sphere to its C++ counterpart using provided properties.
        """
        from PyMieSim.binary.interface_sphere import SPHERE

        self.binding = SPHERE(
            diameter=self.diameter.to_base_units().magnitude,
            refractive_index=self.index.to_base_units().magnitude,
            medium_refractive_index=self.medium_index.to_base_units().magnitude,
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

        Parameters
        ----------
        max_order : Optional[int])
            The maximum order of the coefficient. Default is 0.

        Returns
        -------
        numpy.ndarray
            Array of :math:`a_n` coefficients.
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

        Parameters
        ----------
        max_order : Optional[int])
            The maximum order of the coefficient. Default is 0.

        Returns
        -------
        numpy.ndarray
            Array of :math:`b_n` coefficients.
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

        Parameters
        ----------
        max_order : Optional[int])
            The maximum order of the coefficient. Default is 0.

        Returns
        -------
        numpy.ndarray
            Array of :math:`c_n` coefficients.
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

        Parameters
        ----------
        max_order : Optional[int])
            The maximum order of the coefficient. Default is 0.

        Returns
        -------
        numpy.ndarray
            Array of :math:`d_n` coefficients.
        """
        return self.binding.dn(max_order)
