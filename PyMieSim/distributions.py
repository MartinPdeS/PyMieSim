"""Number-based particle-size distributions for independent-scattering averages."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pint import Quantity

from .units import ureg


def _number_fractions(number_weights: ArrayLike, size: int) -> NDArray[np.float64]:
    """Validate and normalize number counts without overflowing their sum."""
    if isinstance(number_weights, Quantity):
        number_weights = number_weights.to("dimensionless").magnitude
    raw = np.asarray(number_weights)
    if np.iscomplexobj(raw):
        raise ValueError("number_weights must be real")
    weights = np.asarray(raw, dtype=float)
    if weights.shape != (size,) or size == 0:
        raise ValueError("number_weights must be a nonempty one-dimensional array matching the number of sizes or components")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError("number_weights must be finite and nonnegative, with at least one positive weight")
    weights = weights / np.max(weights)
    return weights / np.sum(weights)


def _positive_length(value: Quantity, name: str) -> float:
    """Validate a finite positive scalar length and return meters."""
    if not isinstance(value, Quantity):
        raise TypeError(f"{name} must carry length units")
    magnitude = np.asarray(value.to("meter").magnitude)
    if magnitude.ndim != 0 or np.iscomplexobj(magnitude) or not np.isfinite(magnitude) or magnitude <= 0:
        raise ValueError(f"{name} must be a finite positive scalar")
    return float(magnitude)


def _validate_sampling(sampling: int) -> None:
    """Require a positive number of quadrature points."""
    if isinstance(sampling, bool) or not isinstance(sampling, Integral) or sampling < 1:
        raise ValueError("sampling must be a positive integer")


def _bounds(minimum_diameter: Quantity, maximum_diameter: Quantity) -> tuple[float, float]:
    """Validate ordered, positive diameter bounds."""
    lower = _positive_length(minimum_diameter, "minimum_diameter")
    upper = _positive_length(maximum_diameter, "maximum_diameter")
    if lower >= upper:
        raise ValueError("minimum_diameter must be smaller than maximum_diameter")
    return lower, upper


def _legendre(sampling: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return quadrature nodes and weights on the unit interval."""
    _validate_sampling(sampling)
    nodes, weights = np.polynomial.legendre.leggauss(sampling)
    return (nodes + 1) / 2, weights / 2


@dataclass(frozen=True, init=False)
class ParticleSizeDistribution:
    """Discrete diameters and normalized particle number fractions.

    Parameters
    ----------
    diameters
        One-dimensional, finite, positive diameters with length units.
    number_weights
        Nonnegative particle counts or probability masses at each diameter.
        These are number weights, not volume fractions, intensity weights,
        or samples of a probability density. They are normalized internally.

    Notes
    -----
    This object describes sizes only. ``Experiment.average_size_distribution``
    uses it exclusively in the non-interacting, independent-scattering
    approximation. It does not model interactions between particles, multiple
    scattering, or interparticle interference. For a better approximation when
    particle correlations matter, refer to PackLab's correlation-based
    dependent-scattering calculations:
    https://martinpdes.github.io/PackLab/docs/latest/scattering.html.
    """

    _diameters_m: tuple[float, ...]
    _number_fractions: tuple[float, ...]

    def __init__(self, diameters: Quantity, number_weights: ArrayLike) -> None:
        if not isinstance(diameters, Quantity):
            raise TypeError("diameters must carry length units, for example [100, 200] * ureg.nanometer")
        raw_diameters = np.asarray(diameters.to("meter").magnitude)
        if np.iscomplexobj(raw_diameters):
            raise ValueError("diameters must be real")
        sizes = np.asarray(raw_diameters, dtype=float)
        if sizes.ndim != 1 or sizes.size == 0:
            raise ValueError("diameters must be a nonempty one-dimensional array")
        if not np.all(np.isfinite(sizes)) or np.any(sizes <= 0):
            raise ValueError("diameters must be finite and strictly positive")
        weights = _number_fractions(number_weights, sizes.size)
        object.__setattr__(self, "_diameters_m", tuple(sizes))
        object.__setattr__(self, "_number_fractions", tuple(weights))

    @property
    def diameters(self) -> Quantity:
        """Diameter nodes in meters, returned as a copy."""
        return np.array(self._diameters_m) * ureg.meter

    @property
    def number_fractions(self) -> NDArray[np.float64]:
        """Normalized particle number fractions, returned as a copy."""
        return np.array(self._number_fractions)

    @classmethod
    def lognormal(
        cls, median_diameter: Quantity, geometric_std: float, sampling: int = 32,
    ) -> ParticleSizeDistribution:
        """Discretize a number-based lognormal distribution by Gauss-Hermite quadrature.

        Parameters
        ----------
        median_diameter
            Finite positive scalar diameter with length units. This is the
            number median, not the arithmetic mean or a volume median.
        geometric_std
            Geometric standard deviation, at least one. The standard deviation
            of log(diameter) is ``log(geometric_std)``. One gives a single size.
        sampling
            Number of quadrature nodes (default 32). Increase this and compare
            optical results to verify convergence, especially near resonances.

        Notes
        -----
        The nodes integrate the full lognormal distribution without explicit
        tail truncation. The returned weights are quadrature probability masses,
        not a histogram or pointwise probability density.
        """
        _validate_sampling(sampling)
        if not np.isfinite(geometric_std) or geometric_std < 1:
            raise ValueError("geometric_std must be finite and at least one")
        median = _positive_length(median_diameter, "median_diameter")
        if geometric_std == 1:
            return cls(np.array([float(median)]) * ureg.meter, [1.0])
        nodes, weights = np.polynomial.hermite.hermgauss(sampling)
        with np.errstate(over="ignore", under="ignore"):
            diameters = float(median) * np.exp(np.sqrt(2) * np.log(geometric_std) * nodes)
        return cls(diameters * ureg.meter, weights)

    @classmethod
    def monodisperse(cls, diameter: Quantity) -> ParticleSizeDistribution:
        """Put all particles at one positive diameter with length units."""
        value = _positive_length(diameter, "diameter")
        return cls(np.array([value]) * ureg.meter, [1.0])

    @classmethod
    def uniform(
        cls, minimum_diameter: Quantity, maximum_diameter: Quantity, sampling: int = 32,
    ) -> ParticleSizeDistribution:
        """Discretize a number density uniform in diameter between positive bounds.

        Parameters
        ----------
        minimum_diameter, maximum_diameter
            Finite positive scalar lengths, with minimum strictly below maximum.
        sampling
            Number of Gauss-Legendre quadrature points. Increase it to check
            convergence of optical averages.

        Notes
        -----
        The density is constant in diameter, not in log diameter or volume.
        The returned number fractions include the quadrature weights.
        """
        lower, upper = _bounds(minimum_diameter, maximum_diameter)
        nodes, weights = _legendre(sampling)
        return cls((lower + (upper - lower) * nodes) * ureg.meter, weights)

    @classmethod
    def truncated_normal(
        cls, mean_diameter: Quantity, standard_deviation: Quantity, *,
        minimum_diameter: Quantity, maximum_diameter: Quantity, sampling: int = 64,
    ) -> ParticleSizeDistribution:
        """Discretize a normal number distribution conditioned on positive bounds.

        Parameters
        ----------
        mean_diameter
            Positive scalar mean of the underlying untruncated normal, with
            length units. Truncation generally changes the actual mean.
        standard_deviation
            Positive scalar standard deviation of that underlying normal,
            with length units. Use ``monodisperse`` for zero width.
        minimum_diameter, maximum_diameter
            Required positive scalar bounds. Probability outside this interval
            is discarded and the retained distribution is renormalized.
        sampling
            Number of Gauss-Legendre points, default 64. Increase it to check
            convergence, especially for intervals much wider than the width.

        Notes
        -----
        Weights integrate exp(-0.5*((d-mean)/standard_deviation)**2) on the
        specified interval. They are number fractions, not density samples.
        """
        mean = _positive_length(mean_diameter, "mean_diameter")
        width = _positive_length(standard_deviation, "standard_deviation")
        lower, upper = _bounds(minimum_diameter, maximum_diameter)
        nodes, weights = _legendre(sampling)
        diameters = lower + (upper - lower) * nodes
        with np.errstate(over="ignore", invalid="ignore"):
            log_density = -0.5 * ((diameters - mean) / width) ** 2
        if not np.any(np.isfinite(log_density)):
            raise ValueError("Normal density cannot be resolved on these bounds; check the mean and standard_deviation")
        weights *= np.exp(log_density - np.max(log_density))
        return cls(diameters * ureg.meter, weights)

    @classmethod
    def triangular(
        cls, minimum_diameter: Quantity, mode_diameter: Quantity,
        maximum_diameter: Quantity, sampling: int = 16,
    ) -> ParticleSizeDistribution:
        """Discretize a triangular number density with a specified peak diameter.

        Parameters
        ----------
        minimum_diameter, maximum_diameter
            Finite positive scalar bounds with length units.
        mode_diameter
            Positive scalar peak position, within the closed diameter interval.
            An endpoint mode gives a monotonically increasing or decreasing density.
        sampling
            Gauss-Legendre points per nonempty side of the peak (default 16).
            The result has up to twice this number of nodes. Increase it to
            check convergence of optical averages.

        Notes
        -----
        Integrating each side separately resolves the change of slope at the
        peak. Weights represent particle number, not particle volume.
        """
        lower, upper = _bounds(minimum_diameter, maximum_diameter)
        mode = _positive_length(mode_diameter, "mode_diameter")
        if not lower <= mode <= upper:
            raise ValueError("mode_diameter must lie between the minimum and maximum diameters")
        nodes, weights = _legendre(sampling)
        left_fraction = (mode - lower) / (upper - lower)
        diameters = []
        fractions = []
        if mode > lower:
            diameters.append(lower + (mode - lower) * nodes)
            fractions.append(2 * weights * nodes * left_fraction)
        if mode < upper:
            diameters.append(mode + (upper - mode) * nodes)
            fractions.append(2 * weights * (1 - nodes) * (1 - left_fraction))
        return cls(np.concatenate(diameters) * ureg.meter, np.concatenate(fractions))

    @classmethod
    def mixture(
        cls, distributions: Sequence[ParticleSizeDistribution], number_weights: ArrayLike,
    ) -> ParticleSizeDistribution:
        """Combine distributions with component particle counts or number fractions.

        Parameters
        ----------
        distributions
            Nonempty sequence of size distributions. Components can use
            different families and numbers of quadrature nodes.
        number_weights
            Nonnegative number weights, one per component, normalized internally.
            A component's weight is its fraction of particles, independent of
            the number of nodes used to discretize it.

        Notes
        -----
        Nodes are concatenated in component order; duplicate diameters are
        retained. This combines sizes only, not size-dependent compositions.
        Optical averaging remains in the non-interacting approximation; see
        PackLab for correlation-based dependent-scattering approximations.
        """
        components = tuple(distributions)
        if not components or any(not isinstance(item, cls) for item in components):
            raise ValueError("distributions must be a nonempty sequence of ParticleSizeDistribution objects")
        component_weights = _number_fractions(number_weights, len(components))
        diameters = np.concatenate([item.diameters.magnitude for item in components])
        fractions = np.concatenate([weight * item.number_fractions for weight, item in zip(component_weights, components)])
        return cls(diameters * ureg.meter, fractions)

    def __repr__(self) -> str:
        return f"ParticleSizeDistribution(samples={len(self._diameters_m)}, weighting='number')"


__all__ = ["ParticleSizeDistribution"]
