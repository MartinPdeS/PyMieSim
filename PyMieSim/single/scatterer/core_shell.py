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

    core_diameter: Quantity
    shell_width: Quantity
    source: BaseSource
    core_index: Optional[Quantity] = None
    shell_index: Optional[Quantity] = None
    medium_index: Optional[Quantity] = None
    core_material: Optional[BaseMaterial] = None
    shell_material: Optional[BaseMaterial] = None
    medium_material: Optional[BaseMaterial] = None

    property_names = [
        "size_parameter", "area", "g",
        "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qpr",
        "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cpr"
    ]

    @field_validator('core_diameter', 'shell_width', mode='before')
    def validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units.
        """
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return value

    @field_validator('core_index', 'shell_index', 'medium_index', mode='before')
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
            shell_index=self.shell_index.to_base_units().magnitude,
            core_index=self.core_index.to_base_units().magnitude,
            shell_width=self.shell_width.to_base_units().magnitude,
            core_diameter=self.core_diameter.to_base_units().magnitude,
            medium_index=self.medium_index.to_base_units().magnitude,
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
