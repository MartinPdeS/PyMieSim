

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import RefractiveIndex
from PyMieSim.binary.interface_experiment import ScattererProperties, MediumProperties


def get_refractive_index_instance(
    name: str,
    source,
    refractive_index: RefractiveIndex | BaseMaterial,
):

    if not isinstance(refractive_index, (list, tuple, numpy.ndarray)):
        refractive_index = [refractive_index]

    cpp_class = MediumProperties if name == "medium" else ScattererProperties

    if all(isinstance(item, RefractiveIndex) for item in refractive_index):
        return cpp_class([item.magnitude for item in refractive_index])

    if all(isinstance(item, BaseMaterial) for item in refractive_index):
        values = [
            material.compute_refractive_index(source.wavelength)
            for material in refractive_index
        ]
        return cpp_class(numpy.asarray(values).tolist())

    raise TypeError(
        "refractive_index must contain either RefractiveIndex or BaseMaterial objects."
    )
