#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
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
    scatterer_set: Sphere | Cylinder | CoreShell
    """ Scatterer set instance which defined the ranging paremters to be measured """
    source_set: Gaussian | PlaneWave
    """ Source set instance which defined the ranging paremters to be measured, source is PlaneWave """
    detector_set: Photodiode | LPMode = None
    """ Detector set instance which defined the ranging paremters to be measured, default is None as all parameter do not need a detector """

    def __post_init__(self):
        self.scatterer_set.source_set = self.source_set

        self.scatterer_set.initialize_binding()

        self.binding = CppExperiment()

        self.bind_sets_to_experiment()

        self.generate_x_table()

    def generate_x_table(self) -> None:
        self.x_table = self.source_set.append_to_table(table=[])

        self.x_table = self.scatterer_set.append_to_table(table=self.x_table)

        if self.detector_set is not None:
            self.x_table = self.detector_set.append_to_table(table=self.x_table)

    def bind_sets_to_experiment(self) -> None:
        """
        As the name suggest, methods binds all sets cpp_binding to respective class.

        :returns:   No return
        :rtype:     None
        """
        self.source_set.bind_to_experiment(experiment=self)
        self.scatterer_set.bind_to_experiment(experiment=self)

        if self.detector_set:
            self.detector_set.bind_to_experiment(experiment=self)

    def get(self, measure: Table.XParameter, export_as_numpy: bool = False) -> numpy.ndarray | DataVisual:
        """
        Compute the measure provided and return a DataVisual structured array.

        :param      measure:          The measure
        :type       measure:          Table.XParameter
        :param      export_as_numpy:  If true method returns numpy array
        :type       export_as_numpy:  bool

        :returns:   The data structure.
        :rtype:     DataVisual
        """
        self.y_table = [measure]

        measure_string = f'get_{self.scatterer_set.name}_{measure.name}'

        array = getattr(self.binding, measure_string)()

        if export_as_numpy:
            return array

        x_table = Table.Xtable(self.x_table)

        measure.values = array

        return DataVisual(
            x_table=x_table,
            y=measure
        )


# -
