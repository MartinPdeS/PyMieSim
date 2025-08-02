#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from typing import Union, List, Optional, Tuple
from pydantic.dataclasses import dataclass
from PyMieSim.experiment.detector.base import BaseDetector
from PyMieSim.units import Quantity, nanometer, micrometer, meter, degree, radian
from PyMieSim.experiment.scatterer import Sphere, CoreShell


@dataclass
class NearFieldProbe(BaseDetector):
    """
    Near-field electromagnetic field detector for spherical scatterers.

    This detector computes the electromagnetic fields in the near-field region
    around spherical particles using Mie scattering theory. The implementation
    requires internal field coefficients (cn, dn) which are currently available
    only for spherical and core-shell geometries.

    Parameters
    ----------
    sampling_points : numpy.ndarray
        3D coordinates where fields should be computed, shape (N, 3).
        Units should be in length (nanometer, micrometer, etc.).
    field_components : List[str]
        List of field components to compute. Options:
        ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', '|E|', '|H|', '|E|²', '|H|²']
    coordinate_system : str, optional
        Coordinate system for field computation. Options: 'cartesian', 'spherical'.
        Default is 'cartesian'.
    normalize_fields : bool, optional
        Whether to normalize fields by incident field amplitude. Default is True.

    Attributes
    ----------
    supported_geometries : List[type]
        List of scatterer types that support near-field computation.

    Examples
    --------
    >>> import numpy as np
    >>> from PyMieSim.units import nanometer
    >>>
    >>> # Create sampling grid around particle
    >>> x = np.linspace(-500, 500, 50) * nanometer
    >>> y = np.linspace(-500, 500, 50) * nanometer
    >>> z = np.array([0]) * nanometer
    >>> X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    >>> points = np.stack([X.flatten(), Y.flatten(), Z.flatten()], axis=1)
    >>>
    >>> # Create near-field detector
    >>> detector = NearFieldProbe(
    ...     sampling_points=points,
    ...     field_components=['Ex', 'Ey', 'Ez', '|E|²']
    ... )

    Notes
    -----
    The near-field computation is based on the multipole expansion:

    .. math::
        \\mathbf{E}(\\mathbf{r}) = \\sum_{n=1}^{N_{max}} \\sum_{m=-n}^{n}
        [a_n \\mathbf{M}_{mn}^{(3)} + b_n \\mathbf{N}_{mn}^{(3)}]

    For points inside the sphere, internal coefficients c_n and d_n are used
    with regular spherical Bessel functions.

    References
    ----------
    .. [1] Bohren, C.F. and Huffman, D.R., 1983. Absorption and scattering of
           light by small particles. John Wiley & Sons.
    .. [2] Mishchenko, M.I., Travis, L.D. and Lacis, A.A., 2002. Scattering,
           absorption, and emission of light by small particles.
    """

    sampling_points: np.ndarray
    field_components: List[str]
    coordinate_system: str = 'cartesian'
    normalize_fields: bool = True

    # Class attribute for supported geometries
    supported_geometries = [Sphere, CoreShell]

    def __post_init__(self):
        """Initialize the near-field detector."""
        self._validate_inputs()
        self._generate_binding()
        self._generate_mapping()

    def _validate_inputs(self):
        """Validate input parameters for near-field computation."""
        # Validate sampling points
        if not isinstance(self.sampling_points, np.ndarray):
            raise ValueError("sampling_points must be a numpy array")

        if self.sampling_points.ndim != 2 or self.sampling_points.shape[1] != 3:
            raise ValueError("sampling_points must have shape (N, 3) for 3D coordinates")

        # Validate field components
        valid_components = ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', '|E|', '|H|', '|E|²', '|H|²']
        for component in self.field_components:
            if component not in valid_components:
                raise ValueError(f"Invalid field component '{component}'. "
                               f"Valid options: {valid_components}")

        # Validate coordinate system
        if self.coordinate_system not in ['cartesian', 'spherical']:
            raise ValueError("coordinate_system must be 'cartesian' or 'spherical'")

    def _generate_binding(self):
        """Generate C++ binding for near-field computation."""
        # For now, create placeholder binding
        # This will interface with C++ near-field computation engine
        self.binding = None  # TODO: Implement C++ near-field binding

    def _generate_mapping(self):
        """Generate mapping for data organization."""
        self.mapping = {
            'sampling_points': self.sampling_points,
            'field_components': self.field_components,
            'coordinate_system': self.coordinate_system,
            'normalize_fields': self.normalize_fields
        }

    @classmethod
    def validate_scatterer_compatibility(cls, scatterer) -> bool:
        """
        Check if scatterer is compatible with near-field computation.

        Parameters
        ----------
        scatterer : BaseScatterer
            The scatterer object to validate.

        Returns
        -------
        bool
            True if scatterer supports near-field computation.

        Raises
        ------
        NotImplementedError
            If scatterer type doesn't support near-field computation.
        """
        scatterer_type = type(scatterer)

        if scatterer_type not in cls.supported_geometries:
            supported_names = [geom.__name__ for geom in cls.supported_geometries]
            raise NotImplementedError(
                f"Near-field computation not implemented for {scatterer_type.__name__}. "
                f"Currently supported: {supported_names}. "
                f"Cylinder support requires implementation of cn/dn coefficients."
            )

        return True

    @staticmethod
    def create_field_grid(
        x_range: Tuple[Quantity, Quantity],
        y_range: Tuple[Quantity, Quantity],
        z_range: Tuple[Quantity, Quantity],
        resolution: Quantity
    ) -> np.ndarray:
        """
        Create a regular 3D grid for field sampling.

        Parameters
        ----------
        x_range, y_range, z_range : Tuple[Quantity, Quantity]
            (min, max) ranges for each coordinate axis.
        resolution : Quantity
            Grid spacing in length units.

        Returns
        -------
        numpy.ndarray
            Grid points with shape (N, 3).

        Examples
        --------
        >>> from PyMieSim.units import nanometer
        >>> points = NearFieldProbe.create_field_grid(
        ...     x_range=(-500*nanometer, 500*nanometer),
        ...     y_range=(-500*nanometer, 500*nanometer),
        ...     z_range=(0*nanometer, 0*nanometer),
        ...     resolution=50*nanometer
        ... )
        """
        # Convert to base units for computation
        x_vals = np.arange(x_range[0].magnitude, x_range[1].magnitude, resolution.magnitude)
        y_vals = np.arange(y_range[0].magnitude, y_range[1].magnitude, resolution.magnitude)
        z_vals = np.arange(z_range[0].magnitude, z_range[1].magnitude, resolution.magnitude)

        # Create meshgrid
        X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

        # Stack into coordinate array
        points = np.stack([X.flatten(), Y.flatten(), Z.flatten()], axis=1)

        # Add units back
        unit = x_range[0].units
        return points * unit

    @staticmethod
    def create_radial_grid(
        center: Tuple[Quantity, Quantity, Quantity],
        radial_range: Tuple[Quantity, Quantity],
        angular_resolution: Quantity = 10 * degree,
        radial_resolution: Quantity = 50 * nanometer
    ) -> np.ndarray:
        """
        Create a spherical grid centered on the particle.

        Parameters
        ----------
        center : Tuple[Quantity, Quantity, Quantity]
            (x, y, z) center coordinates.
        radial_range : Tuple[Quantity, Quantity]
            (r_min, r_max) radial sampling range.
        angular_resolution : Quantity, optional
            Angular resolution for theta and phi. Default is 10 degrees.
        radial_resolution : Quantity, optional
            Radial resolution. Default is 50 nanometers.

        Returns
        -------
        numpy.ndarray
            Grid points in Cartesian coordinates, shape (N, 3).
        """
        # Create spherical coordinate arrays
        r_vals = np.arange(
            radial_range[0].magnitude,
            radial_range[1].magnitude,
            radial_resolution.magnitude
        )
        theta_vals = np.arange(0, 180, angular_resolution.to(degree).magnitude) * np.pi / 180
        phi_vals = np.arange(0, 360, angular_resolution.to(degree).magnitude) * np.pi / 180

        # Create meshgrid
        R, THETA, PHI = np.meshgrid(r_vals, theta_vals, phi_vals, indexing='ij')

        # Convert to Cartesian
        X = R * np.sin(THETA) * np.cos(PHI) + center[0].magnitude
        Y = R * np.sin(THETA) * np.sin(PHI) + center[1].magnitude
        Z = R * np.cos(THETA) + center[2].magnitude

        # Stack and add units
        points = np.stack([X.flatten(), Y.flatten(), Z.flatten()], axis=1)
        unit = center[0].units
        return points * unit
