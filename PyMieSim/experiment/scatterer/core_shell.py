#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from typing import List
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import Length, RefractiveIndex

from PyMieSim.binary.interface_experiment import CoreShellSet
from PyMieSim.experiment.scatterer.base import BaseScatterer
from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import Sequential



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
    core_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
        Refractive index or indices of the core.
    shell_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
        Refractive index or indices of the shell.
    medium_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
        BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

    """

    source: BaseSource
    core_diameter: Length
    shell_thickness: Length
    core_refractive_index: List[BaseMaterial] | List[RefractiveIndex]
    shell_refractive_index: List[BaseMaterial] | List[RefractiveIndex]
    medium_refractive_index: List[BaseMaterial] | List[RefractiveIndex]

    available_measure_list = [
        "Qsca",
        "Qext",
        "Qabs",
        "Qratio",
        "Qforward",
        "Qback",
        "Qpr",
        "Csca",
        "Cext",
        "Cabs",
        "Cratio",
        "Cforward",
        "Cback",
        "Cpr",
        "a1",
        "a2",
        "a3",
        "b1",
        "b2",
        "b3",
        "g",
        "coupling",
    ]

    attributes = [
        "core_diameter",
        "shell_thickness",
        "core_refractive_index",
        "shell_refractive_index",
        "medium_refractive_index",
    ]

    def __init__(
        self,
        source: BaseSource,
        core_diameter: Length,
        shell_thickness: Length,
        core_refractive_index: List[BaseMaterial] | List[RefractiveIndex],
        shell_refractive_index: List[BaseMaterial] | List[RefractiveIndex],
        medium_refractive_index: List[BaseMaterial] | List[RefractiveIndex],
    ):
        self.mapping = {}
        self.source = source
        self.core_diameter = np.atleast_1d(core_diameter)
        self.shell_thickness = np.atleast_1d(shell_thickness)
        self.core_refractive_index = np.atleast_1d(core_refractive_index)
        self.shell_refractive_index = np.atleast_1d(shell_refractive_index)
        self.medium_refractive_index = np.atleast_1d(medium_refractive_index)

        mediun_refractive_index = self._add_refractive_index(
            name="medium", refractive_index=self.medium_refractive_index
        )

        core_refractive_index = self._add_refractive_index(
            name="core", refractive_index=self.core_refractive_index
        )

        shell_refractive_index = self._add_refractive_index(
            name="shell", refractive_index=self.shell_refractive_index
        )

        self.binding_kwargs = dict(
            core_diameter=self.core_diameter,
            shell_thickness=self.shell_thickness,
            is_sequential=self.is_sequential,
            medium_refractive_index=mediun_refractive_index,
            core_refractive_index=core_refractive_index,
            shell_refractive_index=shell_refractive_index,
        )

        self.set = CoreShellSet(**self.binding_kwargs)
