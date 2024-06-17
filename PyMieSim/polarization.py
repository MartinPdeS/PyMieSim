#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import Iterable, List, Union, NoReturn
import numpy
from pydantic.dataclasses import dataclass


config_dict = dict(
    kw_only=True,
    slots=True,
    extra='forbid'
)


@dataclass(config=config_dict)
class UnitPolarizationAngle:
    angle: float

    def __post_init__(self) -> NoReturn:
        self.jones_vector = numpy.cos(self.angle * numpy.pi / 180), numpy.sin(self.angle * numpy.pi / 180)
        self.jones_vector = numpy.asarray(self.jones_vector)

    def __repr__(self) -> str:
        return self.angle.__repr__()


class JonesVector:
    """
    Represents a Jones vector for describing the polarization state of light.
    """

    def __init__(self, values: Iterable) -> NoReturn:
        """
        Initialize the Jones vector.

        Parameters:
            - jones_vector: Iterable, a collection that represents the Jones vector components.
        """
        self.values = numpy.atleast_2d(values).astype(complex)

    def __repr__(self) -> str:
        return f"JonesVector({self.values})"

    def __add__(self, other) -> 'JonesVector':
        """
        Add another Jones vector to this Jones vector, combining their polarization states.

        Parameters:
            - other: JonesVector, another Jones vector to add.

        Returns:
            - JonesVector, the combined Jones vector.
        """
        return JonesVector(numpy.vstack((self.values, other.values)))


class RightCircular(JonesVector):
    """
    Represents right circular polarization.
    """

    def __init__(self) -> NoReturn:
        super().__init__([1, 1j])


class LeftCircular(JonesVector):
    """
    Represents left circular polarization.
    """

    def __init__(self) -> NoReturn:
        super().__init__([1, -1j])


@dataclass(config=config_dict)
class Linear(JonesVector):
    """
    Represents linear polarization for a given angle or angles.
    """
    angles: Union[List[float], float]

    def __post_init__(self) -> NoReturn:
        """
        Initialize linear polarization with one or more angles.

        Parameters:
            - angles: float, the angle(s) of polarization in degrees.
        """
        self.angles = numpy.atleast_1d(self.angles)

        self.elements = [
            UnitPolarizationAngle(angle=angle) for angle in self.angles
        ]

        self.jones_vector = numpy.asarray([
            unit.jones_vector for unit in self.elements
        ]).astype(complex)

    def __getitem__(self, idx: int) -> UnitPolarizationAngle:
        return self.elements[idx]

    def __repr__(self):
        return self.elements.__repr__()

# -
