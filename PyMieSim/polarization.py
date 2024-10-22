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
    """
    Base class for polarization types.

    This class serves as a parent class for different types of polarizations such as Jones vectors
    and specific polarization states like right circular and left circular polarizations.

    It provides a foundation for shared functionality across all polarization types.
    """
    pass


@dataclass(config=config_dict, unsafe_hash=True)
class JonesVector(BasePolarization):
    """
    Represents a general Jones vector for polarization.

    A Jones vector is a two-dimensional complex vector used to describe the state of polarization
    of electromagnetic waves. This class handles general forms of linear and circular polarization.

    Parameters
    ----------
    element : list
        A list representing the Jones vector. The list will be converted to a 2D NumPy array
        of complex numbers upon initialization.

    Attributes
    ----------
    element : numpy.ndarray
        A 2D NumPy array of complex numbers representing the Jones vector.

    Methods
    -------
    __post_init__():
        Converts `element` to a 2D NumPy array with complex data type after initialization.

    __iter__():
        Returns an iterator over polarization angles if available, otherwise over the Jones vector elements.

    __add__(other: BasePolarization) -> BasePolarization:
        Combines two Jones vectors (or derived polarizations) by stacking their elements.

    __repr__() -> str:
        Returns a string representation of the polarization vector or angle.

    __str__() -> str:
        Converts the polarization vector or angle to a string format for easy display.
    """
    element: list

    def __post_init__(self):
        self.element = np.atleast_2d(self.element).astype(complex)

    def __iter__(self):
        if hasattr(self, 'angle'):
            return iter(self.angle)
        else:
            return iter(self.element)

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
    Represents right circular polarization in the Jones formalism.

    Right circular polarization is a specific form of circular polarization where the electric
    field vector rotates in a right-handed sense around the direction of propagation.

    Attributes
    ----------
    element : numpy.ndarray
        A Jones vector corresponding to right circular polarization, set to `[1, 1j]`.

    Methods
    -------
    __init__():
        Initializes the Jones vector for right circular polarization.
    """
    def __init__(self) -> None:
        super().__init__(element=[1, 1j])


class LeftCircular(JonesVector):
    """
    Represents left circular polarization in the Jones formalism.

    Left circular polarization is a specific form of circular polarization where the electric
    field vector rotates in a left-handed sense around the direction of propagation.

    Attributes
    ----------
    element : numpy.ndarray
        A Jones vector corresponding to left circular polarization, set to `[1, -1j]`.

    Methods
    -------
    __init__():
        Initializes the Jones vector for left circular polarization.
    """
    def __init__(self) -> None:
        super().__init__(element=[1, -1j])


@dataclass(config=config_dict, unsafe_hash=True)
class Linear(JonesVector):
    """
    Represents linear polarization for a given angle or set of angles.

    Linear polarization occurs when the electric field vector oscillates in a straight line
    along the propagation direction. This class allows for defining linear polarization at any angle
    in the plane perpendicular to the propagation direction.

    Parameters
    ----------
    element : Quantity
        A `Quantity` object representing the angle(s) of linear polarization. This is converted into
        a Jones vector based on the provided angle in degrees.

    Attributes
    ----------
    element : numpy.ndarray
        A 2D NumPy array representing the Jones vector for linear polarization.
    angle : numpy.ndarray
        Array of polarization angles in the appropriate units (degrees).

    Methods
    -------
    __post_init__():
        Converts the input `element` (polarization angle) to a Jones vector representation
        and stores it in `element`. It also stores the polarization angle in degrees.
    """
    element: Quantity

    def __post_init__(self):
        self.angle = np.atleast_1d(self.element.magnitude) * self.element.units
        angle = self.element.to(degree).magnitude

        self.element = np.empty([len(self.angle), 2]).astype(complex)
        self.element[:, 0] = np.cos(np.deg2rad(angle))
        self.element[:, 1] = np.sin(np.deg2rad(angle))
