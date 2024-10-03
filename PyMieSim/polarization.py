import numpy as np
from pydantic.dataclasses import dataclass
from PyMieSim.units import degree, Quantity
from copy import copy

config_dict = dict(
    kw_only=False,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)

class BasePolarization:
    pass

@dataclass(config=config_dict, unsafe_hash=True)
class JonesVector(BasePolarization):
    """
    Represents linear polarization for a given angle or angles.
    """
    element: list

    def __post_init__(self):
        self.element = np.atleast_2d(self.element).astype(complex)

    def __add__(self, other: BasePolarization) -> BasePolarization:
        output = copy(self)

        if hasattr(output, 'angle'):
            del output.angle

        output.element = np.vstack([self.element, other.element])

        if hasattr(other, 'angle') and hasattr(self, 'angle'):
            output.angle = np.hstack([self.angle, other.angle])

        return output

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        if hasattr(self, 'angle'):
            return self.angle.__str__()

        return self.element.__str__()

class RightCircular(JonesVector):
    """
    Represents right circular polarization.
    """
    def __init__(self) -> None:
        super().__init__(element=[1, 1j])

class LeftCircular(JonesVector):
    """
    Represents left circular polarization.
    """
    def __init__(self) -> None:
        super().__init__(element=[1, -1j])

@dataclass(config=config_dict, unsafe_hash=True)
class Linear(JonesVector):
    """
    Represents linear polarization for a given angle or angles.
    """
    element: Quantity

    def __post_init__(self):
        self.angle = np.atleast_1d(self.element.magnitude) * self.element.units
        angle = self.element.to(degree).magnitude

        self.element = np.empty([len(self.angle), 2]).astype(complex)
        self.element[:, 0] = np.cos(np.deg2rad(angle))
        self.element[:, 1] =  np.sin(np.deg2rad(angle))
