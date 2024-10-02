#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pandas as pd
from pydantic.dataclasses import dataclass

from PyMieSim.binary.Experiment import CppExperiment
from PyMieSim.units import AU, meter, watt

from typing import Union, Optional
from PyMieSim.experiment.scatterer import Sphere, Cylinder, CoreShell
from PyMieSim.experiment.detector import Photodiode, CoherentMode
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.utils import plot_dataframe


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

    def get(self, *measure, scale_unit: bool = False, drop_unique_level: bool = True) -> Union[numpy.ndarray]:
        """
        Executes the simulation to compute and retrieve the specified measure.

        Parameters:
            measure (Table): The measure to be computed by the simulation, defined by the user.

        Returns:
            Union[numpy.ndarray]: The computed data in the specified format, either as raw numerical values in a numpy array or structured for visualization.
        """
        measure = set(numpy.atleast_1d(measure))

        if not measure.issubset(self.scatterer.available_measure_list):
            raise ValueError(f"Cannot compute {measure} for {self.scatterer.__class__}")

        is_complex = True if numpy.any([m[0] in ['a', 'b'] for m in measure]) else False

        df = self.generate_dataframe(measure=list(measure), is_complex=is_complex)

        if drop_unique_level:
            df = self.drop_unique_levels(df)

        scatterer_name = self.scatterer.__class__.__name__.lower()
        for m in measure:
            call_str = f'get_{scatterer_name}_{m}'

            array = getattr(self.binding, call_str)()

            match m[0]:
                case 'C':
                    dtype = meter**2
                case 'c':
                    dtype = watt
                case _:
                    dtype = AU


            df = self._set_data(measure=m, dataframe=df, array=array, dtype=dtype, is_complex=is_complex)

        setattr(df.__class__, 'plot_data', plot_dataframe)

        if scale_unit:
            df = self.scale_dataframe_units(df)

        return df

    def _set_data(self, measure: str, dataframe: pd.DataFrame, array: numpy.ndarray, dtype: type, is_complex: bool) -> None:
        """
        Sets the real and imaginary parts of a NumPy array as separate 'real' and 'imag' levels
        in the 'type' index of a pandas DataFrame.

        Parameters:
        ----------
        dataframe : pd.DataFrame
            The target DataFrame with a MultiIndex that includes a 'type' level
            (with possible values 'real' and 'imag').

        array : np.ndarray
            A complex NumPy array whose real and imaginary parts will be assigned to
            the DataFrame.

        dtype : type
            The data type to cast the real and imaginary parts into, using PintArrays.

        is_complex : bool
            A flag that indicates whether the input array is complex. If False, the 'type' level
            is not saved, and the array is treated as purely real.
        """
        if is_complex:
            dataframe[(measure, 'real')] = pd.Series(array.ravel().real, dtype=f'pint[{dtype}]', index=dataframe.index)
            dataframe[(measure, 'imag')] = pd.Series(array.ravel().imag, dtype=f'pint[{dtype}]', index=dataframe.index)

        else:
            dataframe[measure] = pd.Series(array.ravel(), dtype=f'pint[{dtype}]', index=dataframe.index)

        return dataframe

    def scale_dataframe_units(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Scales the units of all columns in a pandas DataFrame to the most compact unit
        based on the maximum value in each column.

        This method iterates over each column in the DataFrame, retrieves the maximum value's
        unit, and converts the entire column to that unit using pint's unit system.

        Parameters:
        ----------
        dataframe : pd.DataFrame
            A DataFrame where each column contains PintArray elements with units.

        Returns:
        -------
        pd.DataFrame
            The updated DataFrame with all columns scaled to the most compact unit
            based on the maximum value in each column.
        """
        max_value_unit = dataframe.max().max().to_compact().units

        for name, col in dataframe.items():
            dataframe[name] = col.pint.to(max_value_unit)

        return dataframe

    def drop_unique_levels(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Drops levels from a MultiIndex DataFrame where only a single unique value is present.

        This method identifies any levels in the DataFrame's index that contain only one
        unique value (i.e., they don't vary across the index) and removes those levels
        to simplify the DataFrame.

        Parameters:
        ----------
        dataframe : pd.DataFrame
            The target DataFrame with a MultiIndex index.

        Returns:
        -------
        pd.DataFrame
            A DataFrame with the unique-value levels removed from the MultiIndex.
        """
        unique_levels = [
            l for l in dataframe.index.names if dataframe.index.get_level_values(l).nunique() <= 1  # must be inferior becauase for whatever reason np.nan from polarization_filter would count to zero!
        ]

        return dataframe.droplevel(unique_levels)

    def generate_dataframe(self, measure, is_complex: bool = False):
        """
        Generates a pandas DataFrame with a MultiIndex based on the mapping
        of 'source' and 'scatterer' from an experiment object.

        Parameters:
        experiment: The experiment object containing 'source' and 'scatterer' mappings.

        Returns:
        pd.DataFrame: A DataFrame with a MultiIndex created from the experiment mappings.
        """
        self._generate_mapping()

        iterables = dict()

        iterables.update(self.source.mapping)
        iterables.update(self.scatterer.mapping)

        if self.detector is not None:
            iterables.update(self.detector.mapping)

        # Create a MultiIndex from the iterables
        row_index = pd.MultiIndex.from_product(
            iterables.values(), names=iterables.keys()
        )

        if is_complex:
            columns = pd.MultiIndex.from_product(
                [measure, ['real', 'imag']],
                names=['data', 'type']
            )
        else:
            columns = pd.MultiIndex.from_product(
                [measure],
                names=['data']
            )

        # Return an empty DataFrame with the generated MultiIndex
        df = pd.DataFrame(columns=columns, index=row_index)

        return df

# -
