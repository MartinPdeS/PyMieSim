"""Python-facing API for single-scatterer simulations."""

from typing import Any

from ..measures import MeasureLike, normalize_measures
from ..results import SimulationResult, SimulationResults
from .setup import Setup


class Simulation:
    """Run a single-scatterer PyMieSim simulation.

    The underlying setup is kept private so the public API remains small.
    """

    def __init__(self, scatterer: Any, source: Any, detector: Any = None, debug_mode: bool = False):
        self._setup = Setup(
            scatterer=scatterer,
            source=source,
            detector=detector,
            debug_mode=debug_mode,
        )

    @property
    def available_measures(self) -> tuple[str, ...]:
        """Measures supported by this simulation configuration."""

        names = tuple(getattr(self._setup.scatterer, "property_names", ()))
        if self._setup.detector is not None and "coupling" not in names:
            names += ("coupling",)
        return names

    def run(self, *measures: MeasureLike, as_result: bool = False, **options: Any):
        """Compute measures using the stable single-simulation interface.

        ``Measure`` members and historical strings are both accepted.
        ``as_result=True`` returns explicit typed result containers.
        """

        names = normalize_measures(measures)
        if not names:
            raise ValueError("At least one measure must be requested.")
        if as_result:
            values = {name: SimulationResult(name, self._setup.get(name, **options)) for name in names}
            return next(iter(values.values())) if len(values) == 1 else SimulationResults(values)

        return self._setup.get(*names, **options)

    def get(self, *measures: MeasureLike, **options: Any):
        """Alias for :meth:`run`."""

        return self.run(*measures, **options)

    def plot_system(
        self,
        show_axes: bool = False,
        show_colorbar: bool = True,
        show_detector_cone: bool = False,
        show_unit_sphere: bool = True,
        figure_size: float = 7.0,
    ):
        """Plot the configured source, scatterer, and detector system."""

        return self._setup.plot_system(
            show_axes=show_axes,
            show_colorbar=show_colorbar,
            show_detector_cone=show_detector_cone,
            show_unit_sphere=show_unit_sphere,
            figure_size=figure_size,
        )

    def __repr__(self) -> str:
        return f"<Simulation setup={self._setup!r}>"


__all__ = ["Simulation"]
