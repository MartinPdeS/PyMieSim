"""Python-facing API for single-scatterer simulations."""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable, Literal, cast, overload

from pint import Quantity
from matplotlib.figure import Figure

from .scatterer import BaseScatterer
from .source import BaseSource
from .detector import BaseDetector

from ..mesh import FullMesh
from ..measures import MeasureLike, validate_measures
from ..results import SimulationResult, SimulationResults
from .setup import Setup

if TYPE_CHECKING:
    from .representations import FarFields, Stokes, SPF, S1S2, NearFields, Footprint


class Simulation:
    """Run a single-scatterer PyMieSim simulation.

    The underlying setup is kept private so the public API remains small.
    """

    def __init__(self, scatterer: BaseScatterer, source: BaseSource, detector: BaseDetector | None = None, debug_mode: bool = False) -> None:
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

    @overload
    def run(self, measure: MeasureLike, /, *, as_result: Literal[True]) -> SimulationResult: ...

    @overload
    def run(self, first: MeasureLike, second: MeasureLike, /, *measures: MeasureLike, as_result: Literal[True]) -> SimulationResults: ...

    @overload
    def run(self, measure: MeasureLike, /, *, as_result: Literal[False] = False) -> Quantity: ...

    @overload
    def run(self, first: MeasureLike, second: MeasureLike, /, *measures: MeasureLike, as_result: Literal[False] = False) -> dict[str, Quantity]: ...

    @overload
    def run(self, *measures: MeasureLike, as_result: bool = False) -> Quantity | dict[str, Quantity] | SimulationResult | SimulationResults: ...

    def run(self, *measures: MeasureLike, as_result: bool = False) -> Quantity | dict[str, Quantity] | SimulationResult | SimulationResults:
        """Compute named measures, accepting strings or ``Measure`` members.

        One measure returns a quantity; multiple measures return an ordered
        dictionary. With ``as_result=True``, both simulations and experiments
        return ``SimulationResult`` or ``SimulationResults``, respectively.
        """
        names = validate_measures(measures, self.available_measures)
        values = {name: self._setup.get(name) for name in names}
        if as_result:
            results = {name: SimulationResult(name, value) for name, value in values.items()}
            return next(iter(results.values())) if len(results) == 1 else SimulationResults(results)
        return next(iter(values.values())) if len(values) == 1 else values

    get = run

    def get_representation(self, representation_type: str, **options: object) -> FarFields | Stokes | SPF | S1S2 | NearFields | Footprint:
        """Build a named single-scatterer field representation."""

        return self._setup.get_representation(representation_type, **options)

    @overload
    def get_farfields(self, sampling: int, distance: Quantity) -> tuple[Quantity, Quantity, FullMesh]: ...

    @overload
    def get_farfields(self, phi: Quantity, theta: Quantity, distance: Quantity) -> tuple[Quantity, Quantity]: ...

    def get_farfields(self, *args: object, **options: object) -> tuple[Quantity, Quantity] | tuple[Quantity, Quantity, FullMesh]:
        """Compute structured or angle-sampled far fields."""

        return cast(Callable[..., tuple[Quantity, Quantity] | tuple[Quantity, Quantity, FullMesh]], self._setup.get_farfields)(*args, **options)

    def get_s1s2(self, angles: Quantity) -> tuple[Quantity, Quantity]:
        """Compute the complex angular scattering amplitudes ``S1`` and ``S2``."""

        return self._setup.get_s1s2(angles=angles)

    @overload
    def get_stokes(self, sampling: int, distance: Quantity) -> tuple[Quantity, Quantity, Quantity, Quantity, FullMesh]: ...

    @overload
    def get_stokes(self, phi: Quantity, theta: Quantity, distance: Quantity) -> tuple[Quantity, Quantity, Quantity, Quantity]: ...

    def get_stokes(self, *args: object, **options: object) -> tuple[Quantity, Quantity, Quantity, Quantity] | tuple[Quantity, Quantity, Quantity, Quantity, FullMesh]:
        """Compute the Stokes parameters at selected angles and distance."""

        return cast(Callable[..., tuple[Quantity, Quantity, Quantity, Quantity] | tuple[Quantity, Quantity, Quantity, Quantity, FullMesh]], self._setup.get_stokes)(*args, **options)

    def plot_system(
        self,
        show_axes: bool = False,
        show_colorbar: bool = True,
        show_detector_cone: bool = False,
        show_unit_sphere: bool = True,
        figure_size: float = 7.0,
    ) -> Figure:
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
