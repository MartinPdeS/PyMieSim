#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from typing import List
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, RefractiveIndex

from PyMieSim.binary.interface_sets import CppCoreShellSet
from PyMieSim.experiment.scatterer.base import BaseScatterer
from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import Sequential
from PyMieSim.utils import config_dict


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
    core_diameter : Length
        Diameters of the core components.
    shell_thickness : Length
        Thicknesses of the shell components.
    core_property : List[BaseMaterial] | List[RefractiveIndex]
        Refractive index or indices of the core.
    shell_property : List[BaseMaterial] | List[RefractiveIndex]
        Refractive index or indices of the shell.
    medium_property : List[BaseMaterial] | List[RefractiveIndex]
        BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

    """
    source: BaseSource
    core_diameter: Length
    shell_thickness: Length
    core_property: List[BaseMaterial] | List[RefractiveIndex]
    shell_property: List[BaseMaterial] | List[RefractiveIndex]
    medium_property: List[BaseMaterial] | List[RefractiveIndex]

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

        mediun_properties = self._add_properties(name='medium', properties=self.medium_property)

        core_properties = self._add_properties(name='core', properties=self.core_property)

        shell_properties = self._add_properties(name='shell', properties=self.shell_property)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.to('meter').magnitude,
            shell_thickness=self.shell_thickness.to('meter').magnitude,
            core_properties=core_properties,
            shell_properties=shell_properties,
            medium_properties=mediun_properties,
            is_sequential=self.is_sequential,
        )
