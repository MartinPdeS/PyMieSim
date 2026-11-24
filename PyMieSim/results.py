"""Unit-aware result containers shared by simulations and parameter sweeps."""

from dataclasses import dataclass, field, replace
from typing import Iterator, Mapping, cast

import numpy as np
from numpy.typing import NDArray
from pint import Quantity, Unit
from matplotlib.axes import Axes
from matplotlib.artist import Artist
from pandas import DataFrame

from .labeled_array import LabeledArray


_MEASURE_LABELS = {
    "Qsca": "Scattering efficiency", "Qext": "Extinction efficiency", "Qabs": "Absorption efficiency",
    "Qback": "Backscattering efficiency", "Qforward": "Forward scattering efficiency",
    "Qpr": "Radiation pressure efficiency", "Qratio": "Backscattering ratio",
    "Csca": "Scattering cross section", "Cext": "Extinction cross section", "Cabs": "Absorption cross section",
    "Cback": "Backscattering cross section", "Cforward": "Forward scattering cross section",
    "Cpr": "Radiation pressure cross section", "g": "Asymmetry parameter", "coupling": "Collected power",
    "cross_section": "Geometric cross section", "size_parameter": "Size parameter",
}


def _axis_summary(*, name: str, axis: NDArray, unit: Unit | None) -> str:
    """Summarize a coordinate range in a common compact unit."""
    if axis.size and axis.dtype.kind in "iuf":
        low, high = float(axis.min()), float(axis.max())
        suffix = ""
        if unit is not None:
            scale = max(abs(low), abs(high))
            compact = (scale * unit).to_compact().units if np.isfinite(scale) and scale else unit
            low, high = (low * unit).to(compact).magnitude, (high * unit).to(compact).magnitude
            suffix = f" {compact:~P}".rstrip()
        span = f"{low:.5g}" if low == high else f"{low:.5g}–{high:.5g}"
        return f"{name}: {span}{suffix} ({axis.size})"
    return f"{name}: {axis.size} values"


@dataclass(frozen=True)
class SimulationResult:
    """One named quantity, optionally carrying experiment dimensions and coordinates.

    Returned by ``Simulation.run`` and ``Experiment.run`` for one measure.
    Unit conversion preserves parameter labels and coordinate units.
    """

    measure: str
    value: Quantity
    dims: tuple[str, ...] = ()
    coords: Mapping[str, NDArray] = field(default_factory=dict)
    coordinate_units: Mapping[str, Unit] = field(default_factory=dict)

    metadata: Mapping[str, object] = field(default_factory=dict)

    def __repr__(self) -> str:
        label = _MEASURE_LABELS.get(self.measure)
        title = f"{label} ({self.measure})" if label else self.measure
        if not self.shape:
            return f"{title} = {self.value:.5g~P}".rstrip()
        axes = [_axis_summary(name=name, axis=np.asarray(self.coords[name]), unit=self.coordinate_units.get(name))
                for name in self.dims if name in self.coords]
        unit_label = f"{self.units:~P}" or "dimensionless"
        return "\n".join([f"{title} [{unit_label}] {self.shape}", *axes])

    @property
    def quantity(self) -> Quantity:
        """Underlying quantity, supporting unit conversion and arithmetic."""
        return self.value

    @property
    def magnitude(self) -> "NDArray | float | complex":
        """Numerical values in the current units."""
        return self.value.magnitude

    @property
    def units(self) -> Unit:
        """Physical units of the result."""
        return cast(Unit, self.value.units)

    def to(self, unit: str | Unit) -> "SimulationResult":
        """Return a converted result with the same parameter coordinates."""
        return replace(self, value=self.value.to(unit))

    def as_labeled_array(self) -> LabeledArray:
        """Return native labeled data for plotting, selection, and tabular export."""
        return LabeledArray(
            np.asarray(self.magnitude), list(self.dims), dict(self.coords),
            {**self.metadata, "units": {self.measure: self.units},
             "coordinate_units": dict(self.coordinate_units),
             "measures": (self.measure,)},
            self.measure,
        )

    @property
    def shape(self) -> tuple[int, ...]:
        """Shape of the numerical values."""
        return np.shape(self.magnitude)

    def as_numpy(self) -> NDArray:
        """Return numerical values in the result's current units."""
        return np.asarray(self.magnitude)

    def as_dataframe(self) -> DataFrame:
        """Return a pandas table with coordinate and measure columns."""
        return self.as_labeled_array().as_dataframe()

    def isel(self, indexers: dict[str, int | slice]) -> "SimulationResult":
        """Select experiment coordinates by integer position or slice."""
        return SimulationResult.from_labeled_array(self.as_labeled_array().isel(indexers))

    def plot(
        self,
        x: str | None = None,
        y: str | None = None,
        ax: Axes | None = None,
        cmap: str = "viridis",
        std: str | None = None,
        alpha: float = 0.4,
        legend: bool = True,
        show: bool = True,
    ) -> list[Artist]:
        """Plot the selected measure against a named parameter dimension.

        Parameters match ``LabeledArray.plot``; ``std`` names a dimension to
        summarize with a shaded standard-deviation band.
        """
        return self.as_labeled_array().plot(
            x=x, y=y, ax=ax, cmap=cmap, std=std, alpha=alpha,
            legend=legend, show=show,
        )

    @classmethod
    def from_labeled_array(cls, data: LabeledArray) -> "SimulationResult":
        """Wrap a single-measure experiment result without losing its grid."""
        if data.name is None or len(data.attrs["measures"]) != 1:
            raise ValueError("Expected a labeled array containing exactly one named measure.")
        return cls(
            data.name, data.as_numpy() * data.attrs["units"][data.name],
            data.dims, dict(data.coords), dict(data.attrs["coordinate_units"]),
            {key: value for key, value in data.attrs.items() if key not in {"units", "coordinate_units", "measures"}},
        )


class SimulationResults(Mapping[str, SimulationResult]):
    """Read-only collection returned for multiple requested typed results."""

    def __init__(self, results: Mapping[str, SimulationResult]) -> None:
        self._results = dict(results)

    def __getitem__(self, key: str) -> SimulationResult:
        return self._results[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._results)

    def __len__(self) -> int:
        return len(self._results)

    def as_labeled_array(self) -> LabeledArray:
        """Combine measures along a named axis for tabular export and plotting.

        Each measure retains its own unit in the array metadata.
        """
        names = tuple(self._results)
        if not names:
            raise ValueError("Cannot combine an empty result collection.")
        first = self._results[names[0]]
        for result in self._results.values():
            if result.dims != first.dims or result.shape != first.shape:
                raise ValueError("All measures must share the same dimensions and shape.")
            for dim in first.dims:
                if (
                    dim not in result.coords
                    or dim not in first.coords
                    or not np.array_equal(result.coords[dim], first.coords[dim])
                    or result.coordinate_units.get(dim) != first.coordinate_units.get(dim)
                ):
                    raise ValueError("All measures must share the same coordinates and coordinate units.")
        return LabeledArray(
            np.stack([result.as_numpy() for result in self._results.values()]),
            ["measure", *first.dims],
            {**dict(first.coords), "measure": np.asarray(names, dtype=object)},
            {**first.metadata, "units": {name: result.units for name, result in self._results.items()},
             "coordinate_units": dict(first.coordinate_units), "measures": names},
            None,
        )

    def as_numpy(self) -> NDArray:
        """Return values stacked in measure order, each in its own units."""
        return self.as_labeled_array().as_numpy()

    def as_dataframe(self) -> DataFrame:
        """Return one table column per measure."""
        return self.as_labeled_array().as_dataframe()

    def __repr__(self) -> str:
        return f"SimulationResults({self._results!r})"


__all__ = ["SimulationResult", "SimulationResults"]
