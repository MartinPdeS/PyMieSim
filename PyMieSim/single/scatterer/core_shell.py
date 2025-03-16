#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.material.base_class import BaseMaterial

import numpy
from typing import Optional
from pydantic.dataclasses import dataclass
from PyMieSim.units import Quantity
from PyMieSim.single.scatterer.base import BaseScatterer, config_dict


@dataclass(config=config_dict, kw_only=True)
class CoreShell(BaseScatterer):
    """
    Class representing a core/shell spherical scatterer.

    Parameters
    ----------
    core_diameter : Quantity
        Diameter of the core of the scatterer.
    shell_thickness : Quantity
        Thickness of the shell surrounding the core.
    core_property : Quantity | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer's core. Only one can be provided.
    shell_property : Quantity | BaseMaterial
        Defines either the refractive index (`Quantity`) or material (`BaseMaterial`) of the scatterer's shell. Only one can be provided.

    """

    core_diameter: Quantity
    shell_thickness: Quantity

    core_property: Quantity | BaseMaterial
    shell_property: Quantity | BaseMaterial

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qpr",
        "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cpr"
    ]

    def __post_init__(self):
        self.core_index, self.core_material = self._assign_index_or_material(self.core_property)

        self.shell_index, self.shell_material = self._assign_index_or_material(self.core_property)

        self.cross_section = numpy.pi * ((self.core_diameter / 2 + self.shell_thickness)) ** 2

        super().__post_init__()

    def set_binding(self) -> None:
        """
        Bind the C++ scatterer class.
        """
        from PyMieSim.binary.interface_coreshell import CORESHELL

        self.binding = CORESHELL(
            shell_refractive_index=self.shell_index.to_base_units().magnitude,
            core_refractive_index=self.core_index.to_base_units().magnitude,
            shell_thickness=self.shell_thickness.to_base_units().magnitude,
            core_diameter=self.core_diameter.to_base_units().magnitude,
            medium_refractive_index=self.medium_index.to_base_units().magnitude,
            source=self.source.binding
        )

    def an(self, max_order: Optional[int] = 0) -> numpy.ndarray:
        r"""
        Compute :math:`a_n` coefficient.

        If max_order is set to zero, the maximum order output is calculated using the Wiscombe criterion.

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
        Compute :math:`b_n` coefficient.

        If max_order is set to zero, the maximum order output is calculated using the Wiscombe criterion.

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
