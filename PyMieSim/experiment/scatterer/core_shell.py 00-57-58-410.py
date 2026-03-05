#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import List
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import Length, RefractiveIndex
from PyMieSim.binary.interface_experiment import CoreShell as CoreShellSet
from PyMieSim.binary.interface_experiment import BaseSourceSet



class CoreShell(CoreShellSet):
    pass
    # def __init__(
    #     self,
    #     source: BaseSourceSet,
    #     core_diameter: Length,
    #     shell_thickness: Length,
    #     core_refractive_index: List[BaseMaterial] | List[RefractiveIndex],
    #     shell_refractive_index: List[BaseMaterial] | List[RefractiveIndex],
    #     medium_refractive_index: List[BaseMaterial] | List[RefractiveIndex],
    # ):
    #     self.source = source

    #     medium_refractive_index = self.get_refractive_index_instance(
    #         name="medium", source=source, refractive_index=medium_refractive_index
    #     )

    #     core_refractive_index = self.get_refractive_index_instance(
    #         name="core", source=source, refractive_index=core_refractive_index
    #     )

    #     shell_refractive_index = self.get_refractive_index_instance(
    #         name="shell", source=source, refractive_index=shell_refractive_index
    #     )


    #     super().__init__(
    #         core_diameter=core_diameter,
    #         shell_thickness=shell_thickness,
    #         medium_refractive_index=medium_refractive_index,
    #         core_refractive_index=core_refractive_index,
    #         shell_refractive_index=shell_refractive_index,
    #         is_sequential=self.is_sequential
    #     )


