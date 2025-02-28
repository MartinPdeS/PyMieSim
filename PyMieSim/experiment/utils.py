import numpy as np
from typing import Union, List, Dict, Optional, TypeVar
from PyMieSim.units import Quantity
from PyOptik.material.base_class import BaseMaterial


T = TypeVar('T')

from pydantic import ConfigDict

# Configuration dictionary for the Pydantic dataclass
config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


def broadcast_params(total_size: Optional[int] = None, **kwargs: Union[T, List[T]]) -> Dict[str, np.ndarray]:
    """
    Broadcast scalar or single-element parameters to NumPy arrays of a given target size.

    Each parameter is first converted to a 1D NumPy array. If a parameter is given as a scalar
    or an array/list of length one, it will be repeated (broadcast) to match the target size.
    The target size is determined by `total_size` if provided, or as the maximum size among all
    parameters otherwise. If any parameter is provided as an array of length greater than one and
    its length does not match the target size, a ValueError is raised.

    Parameters
    ----------
    total_size : int, optional
        The explicit target size for broadcasting. Must be a positive integer if provided.
        If not provided, the target size is determined as the maximum size among the input parameters.
    **kwargs
        Keyword arguments mapping parameter names to values, each of which may be a scalar or a list/array.

    Returns
    -------
    dict
        A dictionary with the same keys as the input parameters, where each value is a NumPy array
        of length equal to the target size.

    Raises
    ------
    ValueError
        If any parameter provided as a list/array has a length greater than one that does not match
        the target size, or if total_size is provided and is less than the maximum size among the inputs.
    AssertionError
        If total_size is provided and is not a positive integer.
    """
    if total_size is not None:
        assert isinstance(total_size, int) and total_size > 0, "total_size must be a positive integer"

    # Convert each parameter to a 1D NumPy array.
    arrays = {name: np.atleast_1d(val) for name, val in kwargs.items()}

    # Determine the target size.
    target_size = total_size or  max(arr.size for arr in arrays.values())

    # Ensure no provided array of size > 1 mismatches the explicit target_size.
    for name, arr in arrays.items():
        if arr.size > 1 and arr.size != target_size:
            raise ValueError(f"Inconsistent sizes: '{name}' has size {arr.size} but expected 1 or {target_size}.")

    # Broadcast any array of size 1 to the target size.
    for name, arr in arrays.items():
        if arr.size == 1:
            arrays[name] = np.full(target_size, arr.item())
            if isinstance(arr, Quantity):
                arrays[name] *= arr.units

    return arrays


class Sequential:
    is_sequential = False

    @classmethod
    def build_sequential(cls, total_size: Optional[int] = None, **kwargs) -> object:
        """
        Broadcast scalar or single-element parameters to NumPy arrays of uniform length
        and construct a new instance of the class.

        Each parameter passed via keyword arguments may be provided as a scalar or as a list/array.
        Scalars or single-element lists/arrays will be repeated to match the length of the longest
        parameter. If any parameter is given as a list/array with more than one element, its length
        must equal the maximum length among all parameters (or the explicit `total_size` if provided);
        otherwise, a ValueError is raised.

        Additionally, if a keyword argument 'source' is provided, it is removed from the
        broadcasted parameters and passed directly to the class constructor.

        Parameters
        ----------
        total_size : int, optional
            The explicit target size for broadcasting. Must be a positive integer if provided.
            If not specified, the target size is determined by the maximum size among the input parameters.
        **kwargs : dict
            Arbitrary keyword arguments mapping parameter names to values. Each value may be a scalar
            or a list/array of values. Examples include parameters such as `source`, `diameter`,
            `property`, or `medium_property`.

        Returns
        -------
        object
            A new instance of the class with all provided parameters broadcasted as NumPy arrays
            of equal length. If 'source' is provided, it is passed directly to the class constructor.

        Raises
        ------
        ValueError
            If any parameter provided as a list/array has a length greater than one that does not
            match the target size.
        """
        is_material = any([
            isinstance(v, BaseMaterial) for v in kwargs.values()
        ])

        assert not is_material, "For the moment no Material can be defined material property for sequential computing, one should use refractive index property. See documentation online."
        source = kwargs.pop("source", None)
        kwargs = broadcast_params(total_size=total_size, **kwargs)

        if source is not None:
            instance = cls(source=source, **kwargs)
        else:
            instance = cls(**kwargs)

        instance.is_sequential = True

        return instance

