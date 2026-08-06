"""Typed, opt-in result containers for the stable Python API."""

from dataclasses import dataclass
from typing import Any, Iterator, Mapping


@dataclass(frozen=True)
class SimulationResult:
    """One named result from a single-scatterer simulation.

    ``value`` is intentionally the original PyMieSim quantity, so unit
    conversion and arithmetic continue to use the project's existing unit
    system.  Use ``as_result=True`` on :meth:`Simulation.run` to opt in.
    """

    measure: str
    value: Any

    @property
    def quantity(self) -> Any:
        """Alias for the underlying unit-aware value."""

        return self.value

    @property
    def magnitude(self) -> Any:
        """Magnitude of the underlying quantity when available."""

        return getattr(self.value, "magnitude", self.value)

    @property
    def units(self) -> Any:
        """Units of the underlying quantity, or ``None`` for plain values."""

        return getattr(self.value, "units", None)

    def to(self, unit: Any) -> "SimulationResult":
        """Return a copy converted to ``unit``."""

        return SimulationResult(self.measure, self.value.to(unit))

    def __repr__(self) -> str:
        return f"SimulationResult(measure={self.measure!r}, value={self.value!r})"


class SimulationResults(Mapping[str, SimulationResult]):
    """Immutable collection returned for a multi-measure typed result."""

    def __init__(self, results: Mapping[str, SimulationResult]):
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
