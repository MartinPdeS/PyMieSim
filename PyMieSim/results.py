"""Unit-aware result containers shared by simulations and parameter sweeps."""

from dataclasses import dataclass, field, replace
from typing import Iterator, Mapping, cast

import numpy as np
from numpy.typing import NDArray
from pint import Quantity, Unit

from .labeled_array import LabeledArray


@dataclass(frozen=True)
class SimulationResult:
    """One named quantity, optionally carrying experiment dimensions and coordinates.

    Use ``as_result=True`` on ``Simulation.run`` or ``Experiment.run``. Unit
    conversion preserves parameter labels and coordinate units.
    """

    measure: str
    value: Quantity
    dims: tuple[str, ...] = ()
    coords: Mapping[str, NDArray] = field(default_factory=dict)
    coordinate_units: Mapping[str, Unit] = field(default_factory=dict)

    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def quantity(self) -> Quantity:
        """Underlying quantity, supporting unit conversion and arithmetic."""
        return self.value

    @property
    def magnitude(self) -> NDArray | float | complex:
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

    def __repr__(self) -> str:
        return f"SimulationResults({self._results!r})"


__all__ = ["SimulationResult", "SimulationResults"]
