#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pandas as pd
from pydantic.dataclasses import dataclass

from DataVisual import Array, Table
from PyMieSim.binary.Experiment import CppExperiment

from typing import Union, Optional
from PyMieSim.experiment.scatterer import Sphere, Cylinder, CoreShell
from PyMieSim.experiment.detector import Photodiode, CoherentMode
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.utils import generate_dataframe, plot_dataframe


@dataclass
class Setup:
    """
    Orchestrates the setup and execution of light scattering experiments using PyMieSim.

    Attributes:
        scatterer (Union[Sphere, Cylinder, CoreShell]): Configuration for the scatterer in the experiment.
            Defines the physical properties of the particle being studied.
        source (Union[Gaussian, PlaneWave]): Configuration for the light source. Specifies the characteristics
            of the light (e.g., wavelength, polarization) illuminating the scatterer.
        detector (Union[Photodiode, CoherentMode, None], optional): Configuration for the detector, if any. Details the
            method of detection for scattered light, including positional and analytical parameters. Defaults to None.

    Methods provide functionality for initializing bindings, generating parameter tables for visualization,
    and executing the simulation to compute and retrieve specified measures.
    """
    scatterer: Union[Sphere, Cylinder, CoreShell]
    source: Union[Gaussian, PlaneWave]
    detector: Optional[Union[Photodiode, CoherentMode]] = None

    def __post_init__(self):
        """
        Initializes the experiment by setting the source for the scatterer and establishing bindings
        between the components and the simulation environment.
        """
        self._initialize_experiment()
        self._bind_components()

    def _initialize_experiment(self) -> None:
        """
        Initializes the experiment with necessary bindings.
        """
        self.scatterer.source = self.source

        self.binding = CppExperiment()

    def _bind_components(self):
        """Binds the experiment components to the CppExperiment instance."""

        self.binding.set_source(self.source.binding)

        method_str = 'set_' + self.scatterer.__class__.__name__.lower()

        getattr(self.binding, method_str)(self.scatterer.binding)

        if self.detector is not None:
            self.binding.set_detector(self.detector.binding)

    def _generate_mapping(self) -> None:
        self.source._generate_mapping()
        self.scatterer._generate_mapping()

        if self.detector:
            self.detector._generate_mapping()

    def get(self, measure: Table, export_as: str = 'dataframe') -> Union[numpy.ndarray, Array]:
        """
        Executes the simulation to compute and retrieve the specified measure.

        Parameters:
            measure (Table): The measure to be computed by the simulation, defined by the user.
            export_as (bool): Determines the format of the returned data. If True, returns a numpy array, otherwise returns a Array object for enhanced visualization capabilities.

        Returns:
            Union[numpy.ndarray, Array]: The computed data in the specified format, either as raw numerical values in a numpy array or structured for visualization with Array.
        """

        if measure.short_label not in self.scatterer.available_measure_list:
            raise ValueError(f"Cannot compute {measure.short_label} for {self.scatterer.__class__.__name__.lower()}")

        measure_string = f'get_{self.scatterer.__class__.__name__.lower()}_{measure.short_label}'

        array = getattr(self.binding, measure_string)()

        self._generate_mapping()

        match export_as.lower():
            case 'dataframe':
                return self._export_as_dataframe(measure, array)
            case 'numpy':
                return self._export_as_numpy(array)

    def _export_as_dataframe(self, measure, array: numpy.ndarray) -> pd.DataFrame:
        df = generate_dataframe(
            experiment=self,
            is_complex=numpy.iscomplexobj(array)
        )

        setattr(df.__class__, 'plot_data', plot_dataframe)

        df[measure.short_label] = array.ravel().view(float)

        return df

    def _export_as_numpy(self, array: numpy.array) -> numpy.array:
        for k, v in self.source.binding_kwargs.items():
            setattr(self.source, k, v)
        for k, v in self.scatterer.binding_kwargs.items():
            setattr(self.scatterer, k, v)
        if self.detector is not None:
            for k, v in self.detector.binding_kwargs.items():
                setattr(self.detector, k, v)

        return array


# -
