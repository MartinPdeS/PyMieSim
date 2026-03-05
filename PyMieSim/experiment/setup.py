#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union, Optional, List, Dict, Iterable
import numpy as np
import pandas as pd

from PyMieSim.units import ureg
from PyMieSim.binary.interface_experiment import SETUP
from PyMieSim.experiment.scatterer import SphereSet, InfiniteCylinderSet, CoreShellSet
from PyMieSim.experiment.detector import PhotodiodeSet, CoherentModeSet
from PyMieSim.experiment.source import GaussianSet, PlaneWaveSet
from PyMieSim.experiment.dataframe_subclass import PyMieSimDataFrame

import PyMieSim


class EmptyDetectorSet(PhotodiodeSet):
    """Sentinel class used when no detector is provided."""
    pass


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

    def __init__(
        self,
        scatterer: Union[SphereSet, InfiniteCylinderSet, CoreShellSet],
        source: Union[GaussianSet, PlaneWaveSet],
        detector: Optional[Union[PhotodiodeSet, CoherentModeSet]] = EmptyDetectorSet(),
    ):
        """
        Parameters
        ----------
        scatterer
            Scatterer parameter set.
        source
            Optical source parameter set.
        detector
            Detector configuration set.
        """

        self.scatterer = scatterer
        self.source = source
        self.detector = detector

        super().__init__(debug_mode=PyMieSim.debug_mode)

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

        method = getattr(self, f"get_{measure}_sequential")

        return method(
            scatterer_set=self.scatterer,
            source_set=self.source,
            detector_set=self.detector,
        )

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

        if "coupling" in measures and isinstance(self.detector, EmptyDetectorSet):
            raise ValueError("Detector must be provided to compute coupling.")

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

            method = getattr(self, f"get_{measure}")

            values = method(
                scatterer_set=self.scatterer,
                source_set=self.source,
                detector_set=self.detector,
            )

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

        mappings.update(self.source.get_mapping())
        mappings.update(self.scatterer.get_mapping())

        if not isinstance(self.detector, EmptyDetectorSet):
            mappings.update(self.detector.get_mapping())

        return mappings

    def _separate_units_and_values(self, mappings: Dict[str, Iterable]):
        """
        Extract numeric values and units from parameter mappings.

        This function also handles PyMieSim C++ Properties objects
        which represent either constant parameters or spectral
        materials.

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
            # PyMieSim Properties objects (C++ bindings)
            # ----------------------------------------------------------
            if hasattr(param_values, "is_constant") and hasattr(param_values, "is_spectral"):

                if param_values.is_constant():
                    values[key] = list(param_values)

                else:
                    # Spectral case
                    # Each material corresponds to one entry in the grid
                    values[key] = list(param_values.material_names)

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
        Construct a DataFrame with a parameter grid defined by the source, scatterer and detector parameter sets.


        Parameters
        ----------
        measures
            Names of the measures to compute.
        drop_unique_level
            Remove parameters that only contain a single value.

        Returns
        -------
        pandas.DataFrame
        """
        mappings = self._collect_parameter_mappings()

        values, units = self._separate_units_and_values(mappings)

        if drop_unique_level:

            values = {k: v for k, v in values.items() if len(v) > 1}

            units = {k: u for k, u in units.items() if k in values}

        parameter_names = list(values.keys())


        mesh = np.meshgrid(*values.values(), indexing="ij")

        flattened = [m.reshape(-1) for m in mesh]

        dataframe = PyMieSimDataFrame(dict(zip(parameter_names, flattened)))

        dataframe.attrs["units"] = units

        for m in measures:
            dataframe[m] = np.nan

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

            method = getattr(self, f"get_{measure}")

            values = method(
                scatterer_set=self.scatterer,
                source_set=self.source,
                detector_set=self.detector,
            )

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
            if isinstance(self.detector, EmptyDetectorSet):
                raise ValueError("Detector required for coupling computation.")
            return ureg.watt

        return ureg.dimensionless