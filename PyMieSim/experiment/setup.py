#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from typing import NoReturn
    from PyMieSim.experiment.scatterer import Sphere, Cylinder, CoreShell
    from PyMieSim.experiment.detector import Photodiode, LPMode
    from PyMieSim.experiment.source import Gaussian, PlaneWave

import numpy
from dataclasses import dataclass

from DataVisual import DataVisual
import DataVisual.tables as Table
from PyMieSim.binary.Experiment import CppExperiment


@dataclass
class Setup(object):
    """
    Orchestrates the setup and execution of light scattering experiments using PyMieSim.

    Attributes:
        scatterer_set (Union[Sphere, Cylinder, CoreShell]): Configuration for the scatterer in the experiment.
            Defines the physical properties of the particle being studied.
        source_set (Union[Gaussian, PlaneWave]): Configuration for the light source. Specifies the characteristics
            of the light (e.g., wavelength, polarization) illuminating the scatterer.
        detector_set (Union[Photodiode, LPMode, None], optional): Configuration for the detector, if any. Details the
            method of detection for scattered light, including positional and analytical parameters. Defaults to None.

    Methods provide functionality for initializing bindings, generating parameter tables for visualization,
    and executing the simulation to compute and retrieve specified measures.
    """
    scatterer_set: Sphere | Cylinder | CoreShell
    source_set: Gaussian | PlaneWave
    detector_set: Photodiode | LPMode | None = None

    def __post_init__(self):
        """
        Initializes the experiment by setting the source for the scatterer and establishing bindings
        between the components and the simulation environment.
        """
        self.initialize_experiment()
        self.bind_components()
        self.generate_x_table()

    def initialize_experiment(self) -> NoReturn:
        """
        Initializes the experiment with necessary bindings.
        """
        self.scatterer_set.source_set = self.source_set
        self.scatterer_set.initialize_binding()
        self.binding = CppExperiment()

    def bind_components(self):
        """Binds the experiment components to the CppExperiment instance."""
        for component in [self.source_set, self.scatterer_set, self.detector_set]:
            if component:
                component.bind_to_experiment(experiment=self)

    def generate_x_table(self) -> NoReturn:
        """
        Generates and populates the 'x_table' with parameters from the source, scatterer, and detector sets.
        This table is instrumental for data visualization and analysis.

        Returns:
            NoReturn
        """
        table = self.source_set.append_to_table(table=[])
        table = self.scatterer_set.append_to_table(table=table)
        self.x_table = table if not self.detector_set else self.detector_set.append_to_table(table=table)

    def get(self, measure: Table.XParameter, export_as_numpy: bool = False) -> numpy.ndarray | DataVisual:
        """
        Executes the simulation to compute and retrieve the specified measure.

        Parameters:
            measure (Table.XParameter): The measure to be computed by the simulation, defined by the user.
            export_as_numpy (bool): Determines the format of the returned data. If True, returns a numpy array,
                                    otherwise returns a DataVisual object for enhanced visualization capabilities.

        Returns:
            Union[numpy.ndarray, DataVisual]: The computed data in the specified format, either as raw numerical
                                              values in a numpy array or structured for visualization with DataVisual.
        """
        measure_string = f'get_{self.scatterer_set.name}_{measure.name}'
        array = getattr(self.binding, measure_string)()

        if export_as_numpy:
            return array

        measure.values = array
        return DataVisual(x_table=Table.Xtable(self.x_table), y=measure)


# -
