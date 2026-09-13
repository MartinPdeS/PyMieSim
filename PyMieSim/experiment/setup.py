#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List, Dict, Iterable, Literal, overload, cast
from pint import Quantity, Unit
import numpy as np

from PyMieSim.units import ureg
from PyMieSim.experiment._setup import Setup as SETUP
from PyMieSim.labeled_array import LabeledArray
from PyMieSim.distributions import ParticleSizeDistribution
from PyMieSim.results import SimulationResult, SimulationResults
from PyMieSim.measures import MeasureLike, Measure, validate_measures



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
        name = validate_measures((measure,), self.available_measures)[0]
        return getattr(self, f"get_{name}_sequential")()

    # ------------------------------------------------------------------
    # Main public interface
    # ------------------------------------------------------------------

    @overload
    def run(self, measure: MeasureLike, /, *, drop_unique_level: bool = True, as_result: Literal[True]) -> SimulationResult: ...

    @overload
    def run(self, first: MeasureLike, second: MeasureLike, /, *measures: MeasureLike, drop_unique_level: bool = True, as_result: Literal[True]) -> SimulationResults: ...

    @overload
    def run(self, *measures: MeasureLike, drop_unique_level: bool = True, as_result: Literal[False] = False) -> LabeledArray: ...

    @overload
    def run(self, *measures: MeasureLike, drop_unique_level: bool = True, as_result: bool = False) -> LabeledArray | SimulationResult | SimulationResults: ...

    def run(
        self, *measures: MeasureLike, drop_unique_level: bool = True, as_result: bool = False,
    ) -> LabeledArray | SimulationResult | SimulationResults:
        """Compute measures while preserving the parameter grid and physical units.

        By default return a native ``LabeledArray``. With ``as_result=True``,
        return the same named result containers as ``Simulation.run``.
        ``drop_unique_level=False`` keeps axes containing a single value.
        """
        names = validate_measures(measures, self.available_measures)
        if as_result:
            results = {
                name: SimulationResult.from_labeled_array(self._build_labeled_array([name], drop_unique_level))
                for name in names
            }
            return next(iter(results.values())) if len(results) == 1 else SimulationResults(results)
        return self._build_labeled_array(names, drop_unique_level)

    get = run

    @overload
    def average_size_distribution(self, distribution: ParticleSizeDistribution, measure: MeasureLike, /, *, drop_unique_level: bool = True, as_result: Literal[True]) -> SimulationResult: ...

    @overload
    def average_size_distribution(self, distribution: ParticleSizeDistribution, first: MeasureLike, second: MeasureLike, /, *measures: MeasureLike, drop_unique_level: bool = True, as_result: Literal[True]) -> SimulationResults: ...

    @overload
    def average_size_distribution(self, distribution: ParticleSizeDistribution, *measures: MeasureLike, drop_unique_level: bool = True, as_result: Literal[False] = False) -> LabeledArray: ...

    @overload
    def average_size_distribution(self, distribution: ParticleSizeDistribution, *measures: MeasureLike, drop_unique_level: bool = True, as_result: bool = False) -> LabeledArray | SimulationResult | SimulationResults: ...

    def average_size_distribution(
        self, distribution: ParticleSizeDistribution, *measures: MeasureLike,
        drop_unique_level: bool = True, as_result: bool = False,
    ) -> LabeledArray | SimulationResult | SimulationResults:
        """Average sphere sizes in the non-interacting approximation only.

        Parameters
        ----------
        distribution
            Number-based distribution with nodes matching the SphereSet diameter
            axis, in the same order. Other parameter axes remain independent.
        measures
            Cross sections (Csca, Cext, Cabs, Cback, Cforward, Cpr), corresponding
            efficiencies, g, or coupling. Amplitudes and ratios are unsupported.
        drop_unique_level
            Remove remaining parameter axes containing only one value.
        as_result
            Return named typed results instead of a LabeledArray. Both retain
            the approximation and weighting rules in their metadata.

        Notes
        -----
        This uses the non-interacting (independent-scattering) approximation:
        particles respond in isolation to the prescribed illumination. It
        excludes interparticle electromagnetic coupling, multiple scattering,
        positional correlations, and interference between different particles.
        For a better approximation when particle correlations matter, refer to
        PackLab's correlation-based dependent-scattering calculations:
        https://martinpdes.github.io/PackLab/docs/latest/scattering.html.
        PackLab is not a general solver for full electromagnetic multiple scattering.

        Cross sections and detector powers are number averages per particle,
        assuming the same illumination for each particle. Efficiencies are
        mean cross section divided by mean projected area; g is weighted by
        scattering cross section. These are not predictions of a dense or
        multiply scattering sample, nor sums of coherent particle fields.
        """
        from .distributions import average_size_distribution

        labeled = average_size_distribution(self, distribution, measures, drop_unique_level)
        if not as_result:
            return labeled
        names = labeled.attrs["measures"]
        if len(names) == 1:
            return SimulationResult.from_labeled_array(labeled)
        results = {}
        for index, name in enumerate(names):
            selected = labeled.isel({"measure": index})
            attrs = {**selected.attrs, "measures": (name,), "units": {name: selected.attrs["units"][name]},
                     "weighting": {name: selected.attrs["weighting"][name]}}
            single = LabeledArray(selected.as_numpy(), list(selected.dims), dict(selected.coords), attrs, name)
            results[name] = SimulationResult.from_labeled_array(single)
        return SimulationResults(results)

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

        result_name = measures[0] if len(measures) == 1 else None
        return LabeledArray(data, dims, coords, attrs, result_name)

    # ------------------------------------------------------------------
    # Measure computation
    # ------------------------------------------------------------------

    def _collect_parameter_mappings(self) -> Dict[str, object]:
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

    def _separate_units_and_values(self, mappings: Dict[str, object]) -> tuple[dict[str, Iterable], dict[str, Unit]]:
        """
        Extract numeric values and units from parameter mappings.

        This function normalizes parameter coordinates for the labeled result,
        including non-numeric simulation objects such as materials.

        Returns
        -------
        Tuple[Dict[str, object], Dict[str, pint.Unit]]
        """

        units: dict[str, Unit] = {}
        values: dict[str, Iterable] = {}

        for key, param_values in mappings.items():

            # ----------------------------------------------------------
            # Pint quantities
            # ----------------------------------------------------------
            if isinstance(param_values, Quantity):
                units[key] = cast(Unit, param_values.units)
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

    def _determine_unit(self, measure: str) -> Unit:
        """Resolve physical units from the shared measure metadata."""
        try:
            kind = Measure(measure).unit_kind
        except ValueError:
            # Native multipole coefficients (a1, b1, ...) are dimensionless.
            kind = "dimensionless"
        return {"area": ureg.meter ** 2, "power": ureg.watt, "dimensionless": ureg.dimensionless}[kind]
