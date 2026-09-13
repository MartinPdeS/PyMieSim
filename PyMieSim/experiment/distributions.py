"""Independent-scattering reduction of spherical-particle parameter grids."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..distributions import ParticleSizeDistribution
from ..labeled_array import LabeledArray
from ..measures import MeasureLike, validate_measures
from .scatterer_set import SphereSet

if TYPE_CHECKING:
    from .setup import Setup


INDEPENDENT_SCATTERING_ASSUMPTIONS = (
    "Non-interacting particles under the same prescribed illumination.",
    "No interparticle electromagnetic coupling or multiple scattering.",
    "No positional correlations or interference between different particles.",
    "Cross sections and coupled powers are averages per particle, not totals for a sample.",
)
CROSS_SECTIONS = ("Csca", "Cext", "Cabs", "Cback", "Cforward", "Cpr")
EFFICIENCIES = {"Q" + name[1:]: name for name in CROSS_SECTIONS}
SUPPORTED_MEASURES = (*CROSS_SECTIONS, *EFFICIENCIES, "g", "coupling")


def average_size_distribution(
    experiment: Setup,
    distribution: ParticleSizeDistribution,
    measures: tuple[MeasureLike, ...],
    drop_unique_level: bool,
) -> LabeledArray:
    """Average validated sphere sweeps using the appropriate physical weights."""
    if not isinstance(distribution, ParticleSizeDistribution):
        raise TypeError("distribution must be a ParticleSizeDistribution")
    if not isinstance(experiment.scatterer_set, SphereSet):
        raise TypeError("Size-distribution averaging currently requires a SphereSet of homogeneous spheres")
    components = (experiment.source_set, experiment.scatterer_set, experiment.detector_set)
    if any(component is not None and component.is_sequential for component in components):
        raise ValueError("Size-distribution averaging requires a Cartesian parameter grid, not sequential sets")
    available = tuple(name for name in SUPPORTED_MEASURES if name in experiment.available_measures)
    names = validate_measures(measures, available)
    dimension = "scatterer:diameter"
    diameters = distribution.diameters.to("meter").magnitude
    configured = experiment.scatterer_set.get_mapping()[dimension].to("meter").magnitude
    if np.shape(configured) != diameters.shape or not np.allclose(configured, diameters, rtol=1e-12, atol=0):
        raise ValueError("SphereSet diameter nodes must match distribution.diameters in the same order")

    cache: dict[str, LabeledArray] = {}
    weights = distribution.number_fractions

    def data(name: str) -> LabeledArray:
        if name not in cache:
            cache[name] = experiment._build_labeled_array([name], drop_unique_level=False)
        return cache[name]

    def mean(name: str) -> NDArray:
        labeled = data(name)
        return np.asarray(np.tensordot(labeled.as_numpy(), weights, axes=(labeled.dims.index(dimension), 0)))

    arrays = []
    weighting = {}
    mean_area = np.sum(weights * (np.pi / 4) * diameters ** 2)
    for name in names:
        if name in EFFICIENCIES:
            arrays.append(mean(EFFICIENCIES[name]) / mean_area)
            weighting[name] = "mean_cross_section / mean_projected_area"
        elif name == "g":
            scattering = data("Csca")
            axis = scattering.dims.index(dimension)
            denominator = mean("Csca")
            if np.any(~np.isfinite(denominator)) or np.any(denominator <= 0):
                raise ValueError("Ensemble g is undefined where the mean scattering cross section is zero or nonpositive")
            numerator = np.tensordot(data("g").as_numpy() * scattering.as_numpy(), weights, axes=(axis, 0))
            arrays.append(numerator / denominator)
            weighting[name] = "scattering_cross_section"
        else:
            arrays.append(mean(name))
            weighting[name] = "number"

    template = next(iter(cache.values()))
    dims = [dim for dim in template.dims if dim != dimension]
    coords = {dim: template.coords[dim] for dim in dims}
    values = arrays[0] if len(arrays) == 1 else np.stack(arrays)
    if len(arrays) > 1:
        dims.insert(0, "measure")
        coords["measure"] = np.asarray(names, dtype=object)
    if drop_unique_level:
        for axis in reversed(range(len(dims))):
            if dims[axis] != "measure" and len(coords[dims[axis]]) == 1:
                values = np.take(values, 0, axis=axis)
                coords.pop(dims.pop(axis))
    attrs = {
        "units": {name: experiment._determine_unit(name) for name in names},
        "coordinate_units": {dim: unit for dim, unit in template.attrs["coordinate_units"].items() if dim in dims},
        "measures": tuple(names),
        "approximation": "independent_scattering",
        "assumptions": INDEPENDENT_SCATTERING_ASSUMPTIONS,
        "weighting": weighting,
        "size_distribution": {"diameters_m": diameters.tolist(), "number_fractions": weights.tolist()},
    }
    return LabeledArray(values, dims, coords, attrs, names[0] if len(names) == 1 else None)
