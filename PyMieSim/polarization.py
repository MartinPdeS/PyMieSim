#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import List, Union, NoReturn, Tuple, Any
import numpy
from pydantic.dataclasses import dataclass
from dataclasses import field


config_dict = dict(
    kw_only=False,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


@dataclass(config=config_dict)
class UnitPolarizationAngle:
    angle: float

    def __post_init__(self) -> NoReturn:
        self.jones_vector = numpy.cos(self.angle * numpy.pi / 180), numpy.sin(self.angle * numpy.pi / 180)
        self.jones_vector = numpy.asarray(self.jones_vector)

    def __repr__(self) -> str:
        return f"{self.angle}"


@dataclass(config=config_dict)
class UnitJonesVector:
    jones_vector: Tuple[Any, Any]

    def __post_init__(self) -> NoReturn:
        self.jones_vector = numpy.asarray(self.jones_vector)

    def __repr__(self) -> str:
        return f"[{self.jones_vector[0]}, {self.jones_vector[1]}]"


@dataclass(config=config_dict)
class BasePolarization:
    elements: List[Union[UnitPolarizationAngle, UnitJonesVector]] = field(init=False)

    @property
    def jones_vector(self) -> numpy.ndarray:
        return numpy.asarray([
            unit.jones_vector for unit in self.elements
        ]).astype(complex)

    def __add__(self, other: object) -> 'BasePolarization':
        """
        Add another Jones vector to this Jones vector, combining their polarization states.

        Parameters:
            - other: BasePolarization, another Jones vector to add.

        Returns:
            - BasePolarization, the combined Jones vector.
        """
        self.elements.extend(other.elements)
        self.__class__ = BasePolarization
        return self

    def __getitem__(self, idx: int) -> Union[UnitPolarizationAngle, UnitJonesVector]:
        return self.elements[idx]


@dataclass(config=config_dict)
class JonesVector(BasePolarization):
    """
    Represents linear polarization for a given angle or angles.
    """
    elements: Union[numpy.ndarray, List[Tuple[Any, Any]], Tuple[Any, Any]]

    def __post_init__(self):
        self.elements = numpy.atleast_2d(self.elements).astype(complex)

        self.elements = [
            UnitJonesVector(element) for element in self.elements
        ]


@dataclass(config=config_dict)
class Linear(BasePolarization):
    """
    Represents linear polarization for a given angle or angles.
    """
    angles: Union[numpy.ndarray, List[float], float]

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


class RightCircular(JonesVector):
    """
    Represents right circular polarization.
    """

    def __init__(self) -> NoReturn:
        super().__init__(elements=[[1, 1j]])


class LeftCircular(JonesVector):
    """
    Represents left circular polarization.
    """

    def __init__(self) -> NoReturn:
        super().__init__(elements=[[1, -1j]])

# -
