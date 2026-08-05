"""Stable, Python-facing entry points for common PyMieSim workflows."""

from __future__ import annotations

from typing import Any, overload

from .measures import Measure, MeasureLike, normalize_measures
from .results import SimulationResult, SimulationResults
from .single import Setup


class Simulation:
    """Run a single-scatterer PyMieSim simulation.

    This facade keeps users independent from generated C++ extension modules.
    The underlying setup remains available through :attr:`setup` for advanced
    workflows.
    """

    def __init__(self, scatterer: Any, source: Any, detector: Any = None, debug_mode: bool = False):
        self._setup = Setup(
            scatterer=scatterer,
            source=source,
            detector=detector,
            debug_mode=debug_mode,
        )

    @property
    def setup(self) -> Setup:
        """Return the legacy backend setup for advanced operations."""

        return self._setup

    @property
    def advanced(self) -> Setup:
        """Explicit alias for the advanced/legacy backend interface."""

        return self._setup

    @property
    def available_measures(self) -> tuple[str, ...]:
        """Measures supported by this simulation configuration."""

        names = tuple(getattr(self._setup.scatterer, "property_names", ()))
        if self._setup.detector is not None and "coupling" not in names:
            names += ("coupling",)
        return names

    @overload
    def run(self, *measures: MeasureLike, as_result: bool = False, **options: Any) -> Any: ...

    def run(self, *measures: MeasureLike, as_result: bool = False, **options: Any):
        """Compute measures using the stable simulation interface.

        ``Measure`` members and historical strings are both accepted.
        ``as_result=True`` opts into explicit typed result containers while
        the default preserves the existing quantity return values.
        """

        names = normalize_measures(measures)
        if not names:
            raise ValueError("At least one measure must be requested.")
        if as_result:
            values = {name: SimulationResult(name, self._setup.get(name, **options)) for name in names}
            return next(iter(values.values())) if len(values) == 1 else SimulationResults(values)

        return self._setup.get(*names, **options)

    def get(self, *measures: MeasureLike, **options: Any):
        """Alias for :meth:`run`.

        Use ``as_result=True`` for explicit typed result containers.
        """

        return self.run(*measures, **options)

    def __getattr__(self, name: str):
        """Preserve legacy access to specialized backend methods."""

        return getattr(self._setup, name)

    def __repr__(self) -> str:
        return f"<Simulation setup={self._setup!r}>"


__all__ = ["Simulation"]
