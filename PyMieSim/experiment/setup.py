#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pandas as pd
from pydantic.dataclasses import dataclass
from PyMieSim.binary.interface_experiment import EXPERIMENT
from PyMieSim.units import AU, meter, watt
from typing import Union, Optional, List
from PyMieSim.experiment.scatterer import Sphere, Cylinder, CoreShell
from PyMieSim.experiment.detector import Photodiode, CoherentMode
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment.dataframe_subclass import PyMieSimDataFrame
import PyMieSim

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
        self.scatterer._generate_binding()

        self.source._generate_binding()

        if self.detector is not None:
            self.detector._generate_binding()

        self.scatterer.source = self.source

        self.binding = EXPERIMENT(debug_mode=PyMieSim.debug_mode)

    def _bind_components(self):
        """Binds the experiment components to the EXPERIMENT instance."""

        self.binding.set_Source(self.source.binding)

        method_str = 'set_' + self.scatterer.__class__.__name__

        getattr(self.binding, method_str)(self.scatterer.binding)

        if self.detector is not None:
            self.binding.set_Detector(self.detector.binding)

    def _generate_mapping(self) -> None:
        self.source._generate_mapping()
        self.scatterer._generate_mapping()

        if self.detector:
            self.detector._generate_mapping()

    def get_sequential(self, measure: str) -> numpy.ndarray:
        """
        Executes the simulation to compute and retrieve the specified measures in a sequential manner contrary to the standard get function.
        This means that ranges are not iterated in a structured fashion into a dataframe but are only run once.
        This methods was developed for specific aims, for the best suited vizualization tools go for the .get() method.

        Parameters
        ----------
        measures : str
            The measure to be computed in the simulation, provided as arguments by the user.
        Returns
        -------
        numpy.ndarray
            An array containing the computed measures.
        """
        scatterer_name = self.scatterer.__class__.__name__
        method_name = f'get_{scatterer_name}_{measure}_sequential'

        # Compute the values using the binding method
        return getattr(self.binding, method_name)()

    def get(self, *measures, scale_unit: bool = False, drop_unique_level: bool = True, add_units: bool = True, as_numpy: bool = False) -> pd.DataFrame:
        """
        Run the simulation to compute specified measures and return the results.

        Parameters
        ----------
        measures : tuple
            Variable-length tuple of measure names to compute, specified by the user.
        scale_unit : bool, optional
            If True, scales the units in the output DataFrame to the most appropriate units. Defaults to False.
        drop_unique_level : bool, optional
            If True, removes levels from the DataFrame's MultiIndex where only one unique value exists, simplifying the DataFrame. Defaults to True.
        add_units : bool, optional
            If True, includes units in the output DataFrame using pint quantities. If False, the output will contain raw numeric values only. Defaults to True.
        as_numpy : bool, optional
            If True, returns the result as a NumPy array instead of a DataFrame. Defaults to False.

        Returns
        -------
        pd.DataFrame or np.ndarray
            A DataFrame containing the computed measures with optional units and simplified MultiIndex, or a NumPy array if `as_numpy` is set to True.
        """
        measures = set(numpy.atleast_1d(measures))

        if 'coupling' in measures:
            assert self.detector is not None, "To compute the coupling power the detector has to be provided to Setup class"

        if as_numpy:
            return self._get_measure_array(measures)

        is_complex = self._is_complex_measure(measures)

        df = self.generate_dataframe(
            measure=list(measures),
            is_complex=is_complex,
            drop_unique_level=drop_unique_level
        )

        df = self._get_measure_dataframe(df, measures, is_complex, add_units)

        # Optionally scale the units in the DataFrame
        if scale_unit:
            df = self.scale_dataframe_units(df)

        return df

    def _is_complex_measure(self, measures: set) -> bool:
        """
        Determines if the measures involve complex values based on measure names.
        No complex value are computed now, as the a & b parameters a return as absolute values for convienience.

        Parameters
        ----------
        measures : set
            A set of measures requested for computation.

        Returns
        -------
        bool
            True if any measure involves complex values, False otherwise.
        """
        return False

    def _get_measure_array(self, measures: List[str]) -> pd.DataFrame:
        """
        Retrieve data for specified measures and return them as a NumPy array.

        Parameters
        ----------
        measures : List[str]
            List of measure names to compute.

        Returns
        -------
        np.ndarray
            A NumPy array containing computed values for the specified measures.
        """
        output_array = []
        for measure in measures:
            scatterer_name = self.scatterer.__class__.__name__
            method_name = f'get_{scatterer_name}_{measure}'

            output_array.append(getattr(self.binding, method_name)())

        return numpy.asarray(output_array).squeeze()

    def _get_measure_dataframe(self, df: pd.DataFrame, measures: List[str], is_complex: bool, add_units: bool) -> pd.DataFrame:
        """
        Populate a DataFrame with computed values for specified measures.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame in which to store computed values for the specified measures.
        measures : List[str]
            List of measure names to compute and add to the DataFrame.
        is_complex : bool
            Indicates whether any of the measures involve complex values.
        add_units : bool
            If True, adds units to the computed values in the DataFrame.

        Returns
        -------
        pd.DataFrame
            Updated DataFrame with computed values for the specified measures, with optional units if `add_units` is True.
        """
        for measure in measures:
            scatterer_name = self.scatterer.__class__.__name__
            method_name = f'get_{scatterer_name}_{measure}'

            # Compute the values using the binding method
            array = getattr(self.binding, method_name)()

            # Determine the unit based on the measure type
            dtype = self._determine_dtype(measure)

            # Set the data in the DataFrame
            df = self._set_data(
                measure=measure,
                dataframe=df,
                array=array,
                dtype=dtype,
                is_complex=is_complex,
                add_units=add_units
            )

        return df

    def _determine_dtype(self, measure: str):
        """
        Determine the appropriate unit (dtype) based on the measure's first character.

        Parameters
        ----------
        measure : str
            The measure type to determine the unit for.

        Returns
        -------
        pint.Quantity
            The appropriate unit for the measure.
        """
        match measure[0]:
            case 'C':
                return meter**2  # Cross-section measure
            case 'c':
                assert self.detector is not None, "Detector needs to be defined in order to measure coupling"
                return watt  # Power measure
            case _:
                return AU  # Arbitrary units for other measures

    def _set_data(self, measure: str, dataframe: pd.DataFrame, array: numpy.ndarray, dtype: type, is_complex: bool, add_units: bool) -> None:
        """
        Sets the real and imaginary parts of a NumPy array as separate 'real' and 'imag' levels
        in the 'type' index of a pandas DataFrame.

        Parameters
        ----------
        dataframe : pd.DataFrame
            The target DataFrame with a MultiIndex that includes a 'type' level (with possible values 'real' and 'imag').
        array : np.ndarray
            A complex NumPy array whose real and imaginary parts will be assigned to the DataFrame.
        dtype : type
            The data type to cast the real and imaginary parts into, using PintArrays.
        is_complex : bool
            A flag that indicates whether the input array is complex. If False, the 'type' level is not saved, and the array is treated as purely real.
        """
        dtype = f'pint[{dtype}]' if add_units else float

        if is_complex:
            dataframe[(measure, 'real')] = pd.Series(array.ravel().real, dtype=dtype, index=dataframe.index)
            dataframe[(measure, 'imag')] = pd.Series(array.ravel().imag, dtype=dtype, index=dataframe.index)

        else:
            dataframe[measure] = pd.Series(array.ravel(), dtype=dtype, index=dataframe.index)

        return dataframe

    def scale_dataframe_units(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Scales the units of all columns in a pandas DataFrame to the most compact unit
        based on the maximum value in each column.

        This method iterates over each column in the DataFrame, retrieves the maximum value's
        unit, and converts the entire column to that unit using pint's unit system.

        Parameters
        ----------
        dataframe : pd.DataFrame
            A DataFrame where each column contains PintArray elements with units.

        Returns
        -------
        pd.DataFrame
            The updated DataFrame with all columns scaled to the most compact unit based on the maximum value in each column.
        """
        max_value_unit = dataframe.max().max().to_compact().units

        for name, col in dataframe.items():
            dataframe[name] = col.pint.to(max_value_unit)

        return dataframe

    def generate_dataframe(self, measure, is_complex: bool = False, drop_unique_level: bool = True):
        """
        Generates a pandas DataFrame with a MultiIndex based on the mapping
        of 'source' and 'scatterer' from an experiment object.

        Parameters
        ----------
        experiment :
            The experiment object containing 'source' and 'scatterer' mappings.

        Returns
        -------
        pd.DataFrame
            A DataFrame with a MultiIndex created from the experiment mappings.
        """
        self._generate_mapping()

        iterables = dict()

        iterables.update(self.source.mapping)

        iterables.update(self.scatterer.mapping)

        if self.detector is not None:
            iterables.update(self.detector.mapping)

        if drop_unique_level:
            _iterables = {
                k: v for k, v in iterables.items() if len(v) > 1
            }

            if len(_iterables) != 0:
                iterables = _iterables

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
            columns = pd.MultiIndex.from_product([measure], names=['data'])

        # Return an empty DataFrame with the generated MultiIndex
        return PyMieSimDataFrame(columns=columns, index=row_index)