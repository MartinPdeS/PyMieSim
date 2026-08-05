"""Stable, Python-facing entry points for common PyMieSim workflows."""

from __future__ import annotations

from typing import Any

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
        """Return the underlying setup for advanced or legacy operations."""

        return self._setup

    def run(self, *measures: str, **options: Any):
        """Compute measures using the standard PyMieSim result interface."""

        return self._setup.get(*measures, **options)

    def get(self, *measures: str, **options: Any):
        """Alias for :meth:`run`."""

        return self.run(*measures, **options)

    def __getattr__(self, name: str):
        """Preserve access to specialized setup methods during the transition."""

        return getattr(self._setup, name)

    def __repr__(self) -> str:
        return f"<Simulation setup={self._setup!r}>"


__all__ = ["Simulation"]
