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
from PyMieSim.labeled_array import LabeledArray
from PyMieSim.measures import Measure, MeasureLike, normalize_measure, normalize_measures



class Setup(SETUP):
    """
    High level orchestration class for PyMieSim experiments.

    This class coordinates the interaction between the source,
    scatterer and detector sets and exposes a simple interface
    for computing simulation measures.

    Structured experiment results are returned as the native C++
    :class:`~PyMieSim.labeled_array.LabeledArray`, preserving parameter
    dimensions, coordinates, and units. Convert explicitly with
    ``result.as_numpy()`` or ``result.as_dataframe()`` when needed.
    """

    # ------------------------------------------------------------------
    # Sequential execution
    # ------------------------------------------------------------------

    @property
    def available_measures(self) -> tuple[str, ...]:
        """Measures supported by the configured scatterer and detector sets."""

        available = list(self.scatterer_set.available_measure_list)
        if self.detector_set is None and "coupling" in available:
            available.remove("coupling")
        return tuple(available)

    def get_sequential(self, measure: MeasureLike) -> np.ndarray:
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
        return getattr(self, f"get_{normalize_measure(measure)}_sequential")()

    # ------------------------------------------------------------------
    # Main public interface
    # ------------------------------------------------------------------

    def get(
        self,
        *measures: MeasureLike,
        drop_unique_level: bool = True,
    ):
        """
        Run the simulation and compute the requested measures.

        Parameters
        ----------
        measures
            Names of the measures to compute.
        drop_unique_level
            Remove parameters that only contain a single value.
        Returns
        -------
        LabeledArray
            Native labeled simulation data with parameter coordinates and units.
        """

        measures = self._normalize_measures(measures)

        return self._build_labeled_array(measures, drop_unique_level)

    def _build_labeled_array(self, measures: List[str], drop_unique_level: bool) -> LabeledArray:
        """Build a native labeled array while preserving the experiment grid."""
        mappings = self._collect_parameter_mappings()
        values, coordinate_units = self._separate_units_and_values(mappings)
        parameter_names = list(values)
        parameter_axes = [np.asarray(values[name]) for name in parameter_names]

        if len(parameter_axes) != len(self.array_shape):
            raise ValueError(
                f"Mismatch between number of parameter axes ({len(parameter_axes)}) "
                f"and setup shape dimensions ({len(self.array_shape)})."
            )

        for name, coordinate, size in zip(parameter_names, parameter_axes, self.array_shape):
            if coordinate.size != size:
                raise ValueError(
                    f"Parameter '{name}' has length {coordinate.size} "
                    f"but corresponding setup axis has size {size}."
                )

        arrays = [np.asarray(getattr(self, f"get_{measure}")()) for measure in measures]
        data = arrays[0] if len(arrays) == 1 else np.stack(arrays, axis=0)

        dims = list(parameter_names)
        coords = {
            name: coordinate
            for name, coordinate in zip(parameter_names, parameter_axes)
        }
        attrs = {
            "units": {measure: self._determine_unit(measure) for measure in measures},
            "coordinate_units": coordinate_units,
            "measures": tuple(measures),
        }

        if len(measures) > 1:
            dims.insert(0, "measure")
            coords["measure"] = np.asarray(measures, dtype=object)

        if drop_unique_level:
            axes_to_drop = [
                axis for axis, coordinate in enumerate(parameter_axes)
                if coordinate.size == 1
            ]
            for axis in reversed(axes_to_drop):
                data = np.take(data, 0, axis=axis + (1 if len(measures) > 1 else 0))
                dims.pop(axis + (1 if len(measures) > 1 else 0))
                coords.pop(parameter_names[axis])

        name = measures[0] if len(measures) == 1 else None
        return LabeledArray(data, dims, coords, attrs, name)

    # ------------------------------------------------------------------
    # Measure computation
    # ------------------------------------------------------------------

    def _normalize_measures(self, measures) -> List[str]:
        """Validate and normalize requested measures while preserving order."""
        normalized = normalize_measures(measures)

        if not normalized:
            raise ValueError("At least one measure must be requested.")

        available = set(self.available_measures)

        invalid = [
            measure
            for measure in normalized
            if not isinstance(measure, str) or measure not in available
        ]
        if invalid:
            available_names = ", ".join(sorted(available))
            raise ValueError(
                f"Unknown measure(s): {invalid}. Available measures: {available_names}."
            )

        return normalized

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

            arrays.append(np.squeeze(np.asarray(values)))

        stacked = np.stack(arrays)
        return stacked[0] if len(arrays) == 1 else stacked

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
