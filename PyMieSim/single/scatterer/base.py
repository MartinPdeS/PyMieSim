#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Tuple
import numpy
from tabulate import tabulate
from TypedUnit import Length, Angle, ureg

from PyMieSim.single import representations


class BaseScatterer:
    """
    A generic class for a scatterer, providing properties and methods to compute various scattering-related quantities.

    """

    def print_properties(self) -> None:
        """Prints a table of the scatterer's properties and their values."""
        data = [getattr(self, name) for name in self.property_names]
        property_dict = {"Property": self.property_names, "Value": data}

        table = tabulate(property_dict, headers="keys")
        print(table)

    def get_farfields_array(
        self, phi: Angle, theta: Angle, r: Length
    ) -> Tuple[numpy.ndarray, numpy.ndarray]:
        """
        Computes the scattering far field for unstructured coordinates.

        The method computes the fields up to a constant phase value.

        Parameters
        ----------
        phi : numpy.ndarray
            The phi angles in radians.
        theta : numpy.ndarray
            The theta angles in radians.
        r : numpy.ndarray
            The radial distances.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]: The computed far fields.
        """
        return self._cpp_get_farfields(
            phi=phi.to("radian").magnitude,
            theta=theta.to("radian").magnitude,
            distance=r.to_base_units().magnitude,
        )

    def get_s1s2_array(self, phi: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
        """Return the S1 and S2 scattering amplitudes for arbitrary ``phi`` angles.

        Parameters
        ----------
        phi : numpy.ndarray
            Array of azimuthal angles in **radians** for which the scattering
            amplitudes are computed.

        Returns
        -------
        Tuple[numpy.ndarray, numpy.ndarray]
            Arrays of ``S1`` and ``S2`` values evaluated at the supplied angles.
        """
        phi = Angle.check(phi)

        return self._cpp_get_s1s2(phi=phi + numpy.pi / 2)

    def get_stokes_array(
        self, phi: Angle, theta: Angle, r: Length = 1 * ureg.meter
    ) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
        """Return the Stokes parameters for arbitrary ``phi`` and ``theta`` angles.

        Parameters
        ----------
        phi : Angle
            Azimuthal angles in radians.
        theta : Angle
            Polar angles in radians. Must be broadcastable with ``phi``.
        r : Length, optional
            Radial distance from the scatterer. Defaults to ``1 * meter``.

        Returns
        -------
        Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]
            The ``I``, ``Q``, ``U`` and ``V`` Stokes parameters.
        """
        phi = Angle.check(phi)
        theta = Angle.check(theta)
        r = Length.check(r)

        phi, theta = numpy.broadcast_arrays(phi, theta)

        E_phi, E_theta = self.get_farfields_array(phi=phi, theta=theta, r=r)

        intensity = numpy.abs(E_phi) ** 2 + numpy.abs(E_theta) ** 2
        I = intensity / numpy.max(intensity)
        Q = (numpy.abs(E_phi) ** 2 - numpy.abs(E_theta) ** 2) / intensity
        U = (+2 * numpy.real(E_phi * E_theta.conjugate())) / intensity
        V = (-2 * numpy.imag(E_phi * E_theta.conjugate())) / intensity

        return I, Q, U, V

    def get_farfield_array(
        self, phi: Angle, theta: Angle, r: Length = 1 * ureg.meter
    ) -> Tuple[numpy.ndarray, numpy.ndarray]:
        """Return the far-field electric fields for arbitrary ``phi`` and ``theta``.

        Parameters
        ----------
        phi : numpy.ndarray
            Azimuthal angles in radians.
        theta : numpy.ndarray
            Polar angles in radians. Must be broadcastable with ``phi``.
        r : Length, optional
            Radial distance from the scatterer. Defaults to ``1 * meter``.

        Returns
        -------
        Tuple[numpy.ndarray, numpy.ndarray]
            Complex ``E_phi`` and ``E_theta`` field components.
        """
        phi = Angle.check(phi)
        theta = Angle.check(theta)
        r = Length.check(r)

        phi, theta = numpy.broadcast_arrays(phi, theta)
        E_phi, E_theta = self.get_farfields_array(phi=phi, theta=theta, r=r)
        return E_phi, E_theta

    def get_nearfield(
        self,
        x_range: tuple[Length, Length] | str = "auto",
        y_range: tuple[Length, Length] | str = "auto",
        z_range: Length | str = "auto",
        sampling: int = 100,
        field_components: list[str] = None,
    ) -> representations.NearField:
        r"""
        Compute near-field electromagnetic fields using internal coefficients cn and dn.

        This method calculates the electromagnetic fields inside and near the scatterer
        using the multipole expansion with vector spherical harmonics. The internal fields
        (r < radius) are computed using cn and dn coefficients, while external fields
        (r > radius) use an and bn coefficients.

        The near-field computation uses the full multipole expansion:

        .. math::
            \vec{E}(\vec{r}) = \sum_{n=1}^{n_{\text{max}}} \left[ c_n \vec{M}_{o1n}^{(1)}(k_1 r) + d_n \vec{N}_{e1n}^{(1)}(k_1 r) \right] \quad \text{(inside)}

        .. math::
            \vec{E}(\vec{r}) = \sum_{n=1}^{n_{\text{max}}} \left[ a_n \vec{M}_{o1n}^{(3)}(k r) + b_n \vec{N}_{e1n}^{(3)}(k r) \right] \quad \text{(outside)}

        Where:

        - :math:`c_n, d_n`: Internal field coefficients (inside the scatterer)
        - :math:`a_n, b_n`: External field coefficients (outside the scatterer)
        - :math:`\vec{M}, \vec{N}`: Vector spherical harmonics
        - :math:`k_1, k`: Wave numbers inside and outside the scatterer

        Parameters
        ----------
        x_range : tuple[Length, Length]
            Range of x coordinates (x_min, x_max) for field computation.
        y_range : tuple[Length, Length]
            Range of y coordinates (y_min, y_max) for field computation.
        z_range : tuple[Length, Length] | Length, optional
            Range of z coordinates (z_min, z_max) for 3D computation, or single z value
            for 2D slice. Default is 0 (xy-plane slice).
        resolution : int | tuple[int, int, int], optional
            Number of points along each axis. If int, same resolution for all axes.
            Default is 100.
        field_components : list[str], optional
            List of field components to compute. Options: ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "|E|", "|H|"].
            Default is all the electric field components.

        Returns
        -------
        representations.NearField
            Near-field representation object with computed field distributions,
            visualization methods, and analysis capabilities.

        Raises
        ------
        RuntimeError
            If cn/dn coefficients are not available for the scatterer type.
        ValueError
            If field components or coordinate ranges are invalid.

        Notes
        -----
        - This method requires that cn and dn coefficients have been computed.
          Currently supports spherical scatterers only. Cylinder support requires
          implementation of cn/dn coefficients for infinite cylinders.
        - For points inside the scatterer (r < radius), uses internal field coefficients.
        - For points outside the scatterer (r >= radius), uses external field coefficients.
        - The computation includes proper vector spherical harmonics and radial functions.

        """
        x_range = [-self.diameter, +self.diameter] if x_range == "auto" else x_range
        y_range = [-self.diameter, +self.diameter] if y_range == "auto" else y_range
        z_range = 0 * ureg.meter if z_range == "auto" else z_range

        # Validate and convert coordinate ranges to base units (meters)
        if isinstance(x_range, (list, tuple)) and len(x_range) == 2:
            x_min = Length.check(x_range[0])
            x_max = Length.check(x_range[1])
            x_range = (x_min, x_max)
        else:
            raise ValueError(
                "x_range must be a tuple of two Quantity values (x_min, x_max)"
            )

        if isinstance(y_range, (list, tuple)) and len(y_range) == 2:
            y_min = Length.check(y_range[0])
            y_max = Length.check(y_range[1])
            y_range = (y_min, y_max)
        else:
            raise ValueError(
                "y_range must be a tuple of two Quantity values (y_min, y_max)"
            )

        # Handle z posiiton
        if isinstance(z_range, Length):
            z_val = Length.check(z_range)
            z_range = z_val
        else:
            raise ValueError("z must be a single Quantity")

        # Set default field components if not specified
        if field_components is None:
            field_components = ["Ex", "Ey", "Ez", "|E|"]

        self._cpp_compute_cn_dn()

        # Create and return NearField representation
        return representations.NearField(
            scatterer=self,
            x_range=x_range,
            y_range=y_range,
            z=z_range,
            sampling=sampling,
            field_components=field_components,
        )


