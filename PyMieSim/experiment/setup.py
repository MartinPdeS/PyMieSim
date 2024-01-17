#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass

from DataVisual import DataV
import DataVisual.tables as Table

from PyMieSim.binary.Experiment import CppExperiment
from PyMieSim.experiment.scatterer import Sphere, Cylinder, CoreShell
from PyMieSim.experiment.detector import Photodiode, LPMode
from PyMieSim.experiment.source import Gaussian


@dataclass
class Setup(object):
    scatterer_set: Sphere | Cylinder | CoreShell
    """ Scatterer set instance which defined the ranging paremters to be measured """
    source_set: Gaussian
    """ Source set instance which defined the ranging paremters to be measured, source is PlaneWave """
    detector_set: Photodiode | LPMode = None
    """ Detector set instance which defined the ranging paremters to be measured, default is None as all parameter do not need a detector """

    def __post_init__(self):
        self.scatterer_set.evaluate_index_material(self.source_set.wavelength)

        self.binding = CppExperiment()

        self.bind_sets_to_experiment()

        self.x_table = self.source_set.append_to_table(table=[])
        self.x_table = self.scatterer_set.append_to_table(table=self.x_table)

        if self.detector_set is not None:
            self.x_table = self.detector_set.append_to_table(table=self.x_table)

    def bind_sets_to_experiment(self):
        self.source_set.bind_to_experiment(self)
        self.scatterer_set.bind_to_experiment(self)

        if self.detector_set:
            self.detector_set.bind_to_experiment(self)

    def Get(self, measure) -> DataV:
        """
        Compute the measure provided and return a DataV structured array.

        :param      measures:  The measures
        :type       measures:  list

        :returns:   The data structure.
        :rtype:     DataV
        """
        self.y_table = [measure]

        measure_string = f'get_{self.scatterer_set.name}_{measure.name}'

        array = getattr(self.binding, measure_string)()

        x_table = Table.Xtable(self.x_table)

        return DataV(
            array=array,
            x_table=x_table,
            y_parameter=measure
        )


# -
