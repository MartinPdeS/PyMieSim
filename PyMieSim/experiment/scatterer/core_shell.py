#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from typing import List
from PyMieSim.binary.interface_sets import CppCoreShellSet
from PyOptik.material.base_class import BaseMaterial
from PyMieSim.units import Quantity
from PyMieSim.experiment.scatterer.base import BaseScatterer
from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import config_dict, Sequential


@dataclass(config=config_dict, kw_only=True)
class CoreShell(BaseScatterer, Sequential):
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
        Diameters of the core components.
    shell_thickness : Quantity
        Thicknesses of the shell components.
    core_property : List[BaseMaterial] | List[Quantity]
        Refractive index or indices of the core.
    shell_property : List[BaseMaterial] | List[Quantity]
        Refractive index or indices of the shell.
    medium_property : List[BaseMaterial] | List[Quantity]
        BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

    """
    source: BaseSource
    core_diameter: Quantity
    shell_thickness: Quantity
    core_property: List[BaseMaterial] | List[Quantity]
    shell_property: List[BaseMaterial] | List[Quantity]
    medium_property: List[BaseMaterial] | List[Quantity]

    available_measure_list = [
        'Qsca', 'Qext', 'Qabs', 'Qratio', 'Qforward', 'Qback', 'Qpr',
        'Csca', 'Cext', 'Cabs', 'Cratio', 'Cforward', 'Cback', 'Cpr', 'a1',
        'a2', 'a3', 'b1', 'b2', 'b3', 'g', 'coupling',
    ]

    attributes = ['core_diameter', 'shell_thickness', 'core_property', 'shell_property', 'medium_property']

    def _generate_binding(self) -> None:
        """
        Assembles the keyword arguments necessary for C++ binding, tailored for core-shell scatterers.
        Prepares structured data from scatterer properties for efficient simulation integration.

        This function populates `binding_kwargs` with values formatted appropriately for the C++ backend used in simulations.
        """
        self.mapping = {}

        self.binding_kwargs = dict(
            core_diameter=self.core_diameter,
            shell_thickness=self.shell_thickness,
            is_sequential=self.is_sequential
        )

        self._add_properties(name='medium', properties=self.medium_property)

        self._add_properties(name='core', properties=self.core_property)

        self._add_properties(name='shell', properties=self.shell_property)

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppCoreShellSet(**binding_kwargs)
