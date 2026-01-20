#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from typing import List
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import Length, RefractiveIndex
from PyMieSim.binary.interface_experiment import CylinderSet
from PyMieSim.experiment.scatterer.base import BaseScatterer
from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import Sequential

class Cylinder(BaseScatterer, Sequential):
    """
    Represents a cylindrical scatterer configuration for PyMieSim simulations.

    Parameters
    ----------
    source : PyMieSim.experiment.source.base.BaseSource
        Light source configuration for the simulation.
    diameter : Length
        Diameter(s) of the cylinder in meters.
    refractive_index : List[BaseMaterial] | List[RefractiveIndex]
        Refractive index or indices of the spherical scatterers themselves.
    medium_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
        BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
    """
    available_measure_list = [
        "Qsca",
        "Qext",
        "Qabs",
        "Csca",
        "Cext",
        "Cabs",
        "a11",
        "a21",
        "a12",
        "a22",
        "a13",
        "a23",
        "b11",
        "b21",
        "b12",
        "b22",
        "b13",
        "b23",
        "coupling",
    ]

    attributes = ["diameter", "refractive_index", "medium_refractive_index"]

    def __init__(
        self,
        source: BaseSource,
        diameter: Length,
        refractive_index: List[BaseMaterial] | List[RefractiveIndex],
        medium_refractive_index: List[BaseMaterial] | List[RefractiveIndex],
    ):
        self.source = source
        self.diameter = np.atleast_1d(diameter)
        self.refractive_index = np.atleast_1d(refractive_index)
        self.medium_refractive_index = np.atleast_1d(medium_refractive_index)

        medium_refractive_index = self._add_refractive_index(
            name="medium", refractive_index=self.medium_refractive_index
        )

        refractive_index = self._add_refractive_index(
            name="scatterer", refractive_index=self.refractive_index
        )

        self.binding_kwargs = dict(
            diameter=self.diameter,
            is_sequential=self.is_sequential,
            medium_refractive_index=medium_refractive_index,
            refractive_index=refractive_index,
        )

        self.set = CylinderSet(**self.binding_kwargs)
