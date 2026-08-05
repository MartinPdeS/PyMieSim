"""Discoverable measure names and metadata for the public API."""

from __future__ import annotations

from enum import Enum
from typing import Final, Iterable, Literal


class Measure(str, Enum):
    """Names of quantities exposed by PyMieSim.

    Enum values intentionally match the historical string API.
    """

    QSCA = "Qsca"
    QEXT = "Qext"
    QABS = "Qabs"
    QBACK = "Qback"
    QFORWARD = "Qforward"
    QRATIO = "Qratio"
    QPR = "Qpr"
    CSCA = "Csca"
    CEXT = "Cext"
    CABS = "Cabs"
    CBACK = "Cback"
    CRATIO = "Cratio"
    CFORWARD = "Cforward"
    CPR = "Cpr"
    CROSS_SECTION = "cross_section"
    SIZE_PARAMETER = "size_parameter"
    G = "g"
    G_WITH_FARFIELDS = "g_with_farfields"
    COUPLING = "coupling"

    @property
    def is_efficiency(self) -> bool:
        """Whether this measure is dimensionless."""

        return self.value.startswith("Q") or self in {
            Measure.G,
            Measure.G_WITH_FARFIELDS,
            Measure.QRATIO,
        }

    @property
    def unit_kind(self) -> str:
        """Stable unit category used by result metadata."""

        if self.value.startswith("C") or self is Measure.CROSS_SECTION:
            return "area"
        if self is Measure.COUPLING:
            return "power"
        return "dimensionless"


MeasureName = Literal[
    "Qsca", "Qext", "Qabs", "Qback", "Qforward", "Qratio", "Qpr",
    "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cforward", "Cpr",
    "cross_section", "size_parameter", "g", "g_with_farfields", "coupling",
]
MeasureLike = str | Measure


Qsca: Final[str] = Measure.QSCA.value
Qext: Final[str] = Measure.QEXT.value
Qabs: Final[str] = Measure.QABS.value
Qback: Final[str] = Measure.QBACK.value
Qforward: Final[str] = Measure.QFORWARD.value
Qratio: Final[str] = Measure.QRATIO.value
Qpr: Final[str] = Measure.QPR.value
Csca: Final[str] = Measure.CSCA.value
Cext: Final[str] = Measure.CEXT.value
Cabs: Final[str] = Measure.CABS.value
Cback: Final[str] = Measure.CBACK.value
Cratio: Final[str] = Measure.CRATIO.value
Cforward: Final[str] = Measure.CFORWARD.value
Cpr: Final[str] = Measure.CPR.value
cross_section: Final[str] = Measure.CROSS_SECTION.value
size_parameter: Final[str] = Measure.SIZE_PARAMETER.value
g: Final[str] = Measure.G.value
g_with_farfields: Final[str] = Measure.G_WITH_FARFIELDS.value
coupling: Final[str] = Measure.COUPLING.value

ALL_MEASURES: Final[tuple[str, ...]] = tuple(measure.value for measure in Measure)


def normalize_measure(measure: MeasureLike) -> str:
    """Return a backend-compatible measure name."""

    return measure.value if isinstance(measure, Measure) else measure


def normalize_measures(measures: Iterable[MeasureLike]) -> list[str]:
    """Normalize an iterable of enum members and strings."""

    return [normalize_measure(measure) for measure in measures]


__all__ = [
    "ALL_MEASURES", "Cabs", "Cback", "Cext", "Cforward", "Cpr", "Cratio",
    "Csca", "Measure", "Qabs", "Qback", "Qext", "Qforward", "Qpr", "Qratio",
    "Qsca", "coupling", "cross_section", "g", "g_with_farfields", "MeasureLike", "MeasureName",
    "normalize_measure", "normalize_measures", "size_parameter",
]
