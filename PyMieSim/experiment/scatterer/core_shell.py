#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from typing import List
from PyMieSim.binary.SetsInterface import CppCoreShellSet
from PyOptik.material.base_class import BaseMaterial
from PyMieSim.units import Quantity
from PyMieSim.experiment.scatterer.base import BaseScatterer, config_dict


@dataclass(config=config_dict, kw_only=True)
class CoreShell(BaseScatterer):
    """
    A data class representing a core-shell scatterer configuration used in PyMieSim simulations.

    This class facilitates the setup and manipulation of core-shell scatterers by providing structured
    attributes and methods that ensure the scatterers are configured correctly for simulations.
    It extends the BaseScatterer class, adding specific attributes and methods relevant to core-shell geometries.

    Parameters
    ----------
    source : PyMieSim.experiment.source.base.BaseSource
        Light source configuration for the simulation.
    core_diameter : Quantity
        Diameters of the core components in meters.
    shell_width : Quantity
        Thicknesses of the shell components in meters.
    core_property : List[BaseMaterial] | List[Quantity]
        Refractive index or indices of the core.
    shell_property : List[BaseMaterial] | List[Quantity]
        Refractive index or indices of the shell.
    medium_property : List[BaseMaterial] | List[Quantity]
        BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

    """
    core_diameter: Quantity
    shell_width: Quantity

    core_property: List[BaseMaterial] | List[Quantity]
    shell_property: List[BaseMaterial] | List[Quantity]

    available_measure_list = [
        'Qsca', 'Qext', 'Qabs', 'Qratio', 'Qforward', 'Qback', 'Qpr',
        'Csca', 'Cext', 'Cabs', 'Cratio', 'Cforward', 'Cback', 'Cpr', 'a1',
        'a2', 'a3', 'b1', 'b2', 'b3', 'g', 'coupling',
    ]

    def __post_init__(self) -> None:
        """
        Assembles the keyword arguments necessary for C++ binding, tailored for core-shell scatterers.
        Prepares structured data from scatterer properties for efficient simulation integration.

        This function populates `binding_kwargs` with values formatted appropriately for the C++ backend used in simulations.
        """
        self.mapping = {}

        self.binding_kwargs = dict(core_diameter=self.core_diameter, shell_width=self.shell_width)

        self._assign_index_or_material(element='core_', property=self.core_property)
        self._assign_index_or_material(element='shell_', property=self.shell_property)
        self._assign_index_or_material(element='medium_', property=self.medium_property)

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppCoreShellSet(**binding_kwargs)
