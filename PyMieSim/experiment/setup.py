#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List, Dict, Iterable
import numpy as np
import pandas as pd

from PyMieSim.units import ureg
from PyMieSim.experiment._setup import Setup as SETUP
from PyMieSim.experiment.scatterer_set import SphereSet, InfiniteCylinderSet, CoreShellSet
from PyMieSim.experiment.detector_set import PhotodiodeSet, CoherentModeSet
from PyMieSim.experiment.source_set import GaussianSet, PlaneWaveSet
from PyMieSim.experiment.dataframe_subclass import PyMieSimDataFrame
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.material_set import MaterialSet
from PyMieSim.material import ConstantMaterial, ConstantMedium



class Setup(SETUP):
    """
    High level orchestration class for PyMieSim experiments.

    This class coordinates the interaction between the source,
    scatterer and detector sets and exposes a simple interface
    for computing simulation measures.

    The class supports two execution modes:

    1. Sequential execution returning raw arrays.
    2. Structured execution returning a pandas DataFrame where
       simulation parameters define the parameter grid.

    Units are stored in ``DataFrame.attrs["units"]`` and the
    DataFrame itself only contains pure numerical values.
    """

    # ------------------------------------------------------------------
    # Sequential execution
    # ------------------------------------------------------------------

    def get_sequential(self, measure: str) -> np.ndarray:
        """
        Compute a measure once using the current parameter sets.

        Unlike :func:`get`, this method does not construct a parameter
        grid and simply returns the raw output of the simulation.

        Parameters
        ----------
        measure
            Name of the measure to compute.

        Returns
        -------
        numpy.ndarray
            Computed values.
        """
        return getattr(self, f"get_{measure}_sequential")()

    # ------------------------------------------------------------------
    # Main public interface
    # ------------------------------------------------------------------

    def get(
        self,
        *measures: str,
        drop_unique_level: bool = True,
        add_units: bool = True,
        as_numpy: bool = False,
        scale_unit: bool = True,
    ):
        """
        Run the simulation and compute the requested measures.

        Parameters
        ----------
        measures
            Names of the measures to compute.
        drop_unique_level
            Remove parameters that only contain a single value.
        add_units
            Store units inside ``DataFrame.attrs["units"]``.
        as_numpy
            Return raw NumPy arrays instead of a DataFrame.
        scale_unit
            Automatically convert results to compact units.

        Returns
        -------
        pandas.DataFrame or numpy.ndarray
        """

        measures = list(set(np.atleast_1d(measures)))

        if as_numpy:
            return self._compute_measure_arrays(measures)

        dataframe = self._build_dataframe(measures, drop_unique_level)

        self._populate_measure_columns(dataframe, measures, add_units)

        if scale_unit:
            dataframe = dataframe.to_compact()

        return dataframe

    # ------------------------------------------------------------------
    # Measure computation
    # ------------------------------------------------------------------

    def _compute_measure_arrays(self, measures: List[str]) -> np.ndarray:
        """
        Return measures as a stacked NumPy array.

        Parameters
        ----------
        measures
            Names of the measures to compute.

        Returns
        -------
        numpy.ndarray
            Computed values.
        """
        arrays = []

        for measure in measures:
            values = getattr(self, f"get_{measure}")()

            arrays.append(np.asarray(values))

        return np.squeeze(np.asarray(arrays))

    # ------------------------------------------------------------------
    # DataFrame generation
    # ------------------------------------------------------------------

    def _collect_parameter_mappings(self) -> Dict[str, Iterable]:
        """
        Collect parameter mappings from source, scatterer and detector.

        Returns
        -------
        dict
            A dictionary containing parameter names as keys and their corresponding values as lists.
        """

        mappings = {}

        mappings.update(self.source_set.get_mapping())
        mappings.update(self.scatterer_set.get_mapping())

        if self.detector_set is not None:
            mappings.update(self.detector_set.get_mapping())

        return mappings

    def _separate_units_and_values(self, mappings: Dict[str, Iterable]):
        """
        Extract numeric values and units from parameter mappings.

        This function also converts non-numeric simulation objects such as
        materials into dataframe-safe representations for grouping and plotting.

        Returns
        -------
        Tuple[Dict[str, Iterable], Dict[str, pint.Unit]]
        """

        units = {}
        values = {}

        for key, param_values in mappings.items():

            # ----------------------------------------------------------
            # Pint quantities
            # ----------------------------------------------------------
            if hasattr(param_values, "units"):
                units[key] = param_values.units
                values[key] = np.atleast_1d(param_values.magnitude)
                continue

            # ----------------------------------------------------------
            # Standard iterable parameters
            # ----------------------------------------------------------
            if isinstance(param_values, (list, tuple, np.ndarray)):
                values[key] = list(param_values)
                continue

            # ----------------------------------------------------------
            # Scalar fallback
            # ----------------------------------------------------------
            values[key] = [param_values]

        return values, units

    def _build_dataframe(self, measures: List[str], drop_unique_level: bool):
        """
        Construct a DataFrame from the canonical experiment shape stored in ``self.shape``.

        Parameters
        ----------
        measures
            Names of the measures to initialize in the dataframe.
        drop_unique_level
            Remove parameters that only contain a single value.

        Returns
        -------
        PyMieSimDataFrame
        """
        mappings = self._collect_parameter_mappings()
        values, units = self._separate_units_and_values(mappings)

        parameter_names = list(values.keys())
        parameter_axes = list(values.values())

        if len(parameter_axes) != len(self.array_shape):
            raise ValueError(
                f"Mismatch between number of parameter axes ({len(parameter_axes)}) "
                f"and setup shape dimensions ({len(self.array_shape)})."
            )

        dataframe_dict = {}

        total_size = int(np.prod(self.array_shape))

        for axis_index, (parameter_name, axis_values, axis_size) in enumerate(zip(parameter_names, parameter_axes, self.array_shape)):

            if len(axis_values) != axis_size:
                raise ValueError(
                    f"Parameter '{parameter_name}' has length {len(axis_values)} but corresponding setup axis has size {axis_size}."
                )

            if drop_unique_level and axis_size == 1:
                continue

            axis_object_array = np.empty(axis_size, dtype=object)
            axis_object_array[:] = list(axis_values)

            reshaped = axis_object_array.reshape(
                [axis_size if i == axis_index else 1 for i in range(len(self.array_shape))]
            )

            broadcasted = np.broadcast_to(reshaped, self.array_shape)

            dataframe_dict[parameter_name] = broadcasted.reshape(total_size)

        if drop_unique_level:
            units = {key: value for key, value in units.items() if key in dataframe_dict}

        dataframe = PyMieSimDataFrame(dataframe_dict)

        dataframe.attrs["units"] = units

        for measure in measures:
            dataframe[measure] = np.nan

        return dataframe

    # ------------------------------------------------------------------
    # Populate simulation outputs
    # ------------------------------------------------------------------

    def _populate_measure_columns(
        self,
        dataframe: pd.DataFrame,
        measures: List[str],
        add_units: bool,
    ):
        """
        Fill the DataFrame with computed simulation results.

        Parameters
        ----------
        dataframe
            DataFrame to populate.
        measures
            List of measures to compute and add to the DataFrame.
        add_units
            Whether to store units in ``DataFrame.attrs["units"]``.

        Returns
        -------
        None
        """

        units = dataframe.attrs.setdefault("units", {})

        for measure in measures:

            values = getattr(self, f"get_{measure}")()

            dataframe[measure] = values.ravel()

            if add_units:
                units[measure] = self._determine_unit(measure)

    # ------------------------------------------------------------------
    # Unit inference
    # ------------------------------------------------------------------

    def _determine_unit(self, measure: str):
        """
        Infer the physical unit associated with a measure.

        Parameters
        ----------
        measure
            Name of the measure.

        Returns
        -------
        pint.Unit

        """
        if measure.startswith("C"):
            return ureg.meter ** 2

        if measure.startswith("c"):
            return ureg.watt

        return ureg.dimensionless