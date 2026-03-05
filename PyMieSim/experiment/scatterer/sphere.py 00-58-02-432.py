#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import List
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import Length, RefractiveIndex
from PyMieSim.binary.interface_experiment import Sphere as SphereSet
from PyMieSim.binary.interface_experiment import BaseSourceSet


class Sphere(SphereSet):
    pass
    # def __init__(
    #     self,
    #     source: BaseSourceSet,
    #     diameter: Length,
    #     refractive_index: List[BaseMaterial] | List[RefractiveIndex],
    #     medium_refractive_index: List[BaseMaterial] | List[RefractiveIndex]
    # ):
    #     medium_refractive_index = self.get_refractive_index_instance(
    #         name="medium", source=source, refractive_index=medium_refractive_index
    #     )

    #     refractive_index = self.get_refractive_index_instance(
    #         name="scatterer", source=source, refractive_index=refractive_index
    #     )

    #     super().__init__(
    #         diameter=diameter,
    #         medium_refractive_index=medium_refractive_index,
    #         refractive_index=refractive_index,
    #         is_sequential=self.is_sequential
    #     )
