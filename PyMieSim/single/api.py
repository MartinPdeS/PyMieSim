"""Python-facing API for single-scatterer simulations."""

from typing import Callable, cast, overload

from pint import Quantity
from matplotlib.figure import Figure

from .scatterer import BaseScatterer
from .source import BaseSource
from .detector import BaseDetector

from ..mesh import FullMesh
from ..measures import MeasureLike, validate_measures
from ..results import SimulationResult, SimulationResults
from .setup import Setup
from .representations import FarFields, Stokes, SPF, S1S2, NearFields, Footprint
from PyMieSim.units import ureg



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
    def run(self, measure: MeasureLike, /) -> SimulationResult: ...

    @overload
    def run(self, first: MeasureLike, second: MeasureLike, /, *measures: MeasureLike) -> SimulationResults: ...

    def run(self, *measures: MeasureLike) -> SimulationResult | SimulationResults:
        """Compute one or more named measures as unit-aware results."""
        names = validate_measures(measures, self.available_measures)
        results = {name: SimulationResult(name, self._setup.get(name)) for name in names}
        return next(iter(results.values())) if len(results) == 1 else SimulationResults(results)

    def farfields(self, *, sampling: int = 200, distance: Quantity = 1.0 * ureg.meter) -> FarFields:
        """Compute the angular far-field representation at the given distance."""
        return FarFields(setup=self._setup, sampling=sampling, distance=distance)

    def stokes(self, *, sampling: int = 200, distance: Quantity = 1.0 * ureg.meter) -> Stokes:
        """Compute Stokes parameters on a spherical angular mesh."""
        return Stokes(setup=self._setup, sampling=sampling, distance=distance)

    def spf(self, *, sampling: int = 200) -> SPF:
        """Compute the scattering phase function on an angular mesh."""
        return SPF(setup=self._setup, sampling=sampling)

    def s1s2(self, *, sampling: int = 200) -> S1S2:
        """Compute the angular scattering amplitude representation."""
        return S1S2(setup=self._setup, sampling=sampling)

    def nearfields(self) -> NearFields:
        """Create a near-field representation for the configured scatterer."""
        return NearFields(setup=self._setup)

    def footprint(self, *, sampling: int = 200, padding_factor: int = 20) -> Footprint:
        """Compute the detector-plane footprint for a configured detector."""
        return Footprint(setup=self._setup, sampling=sampling, padding_factor=padding_factor)

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
