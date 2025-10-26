#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Tuple
import numpy
from tabulate import tabulate
from PyOptik.material.base_class import BaseMaterial
from TypedUnit import RefractiveIndex, Length, Dimensionless, Area, Angle, ureg

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

    @property
    def size_parameter(self) -> Dimensionless:
        """Returns the size parameter of the scatterer."""
        return self._cpp_size_parameter * ureg.AU

    @property
    def cross_section(self) -> Dimensionless:
        """Returns the cross-section of the scatterer."""
        return (self._cpp_cross_section * ureg.meter**2).to_compact()

    @property
    def Qsca(self) -> Dimensionless:
        """Returns the scattering efficiency."""
        return self._cpp_Qsca * ureg.AU

    @property
    def Qext(self) -> Dimensionless:
        """Returns the extinction efficiency."""
        return self._cpp_Qext * ureg.AU

    @property
    def Qabs(self) -> Dimensionless:
        """Returns the absorption efficiency."""
        return self._cpp_Qabs * ureg.AU

    @property
    def Qback(self) -> Dimensionless:
        """Returns the backscattering efficiency."""
        return self._cpp_Qback * ureg.AU

    @property
    def Qforward(self) -> Dimensionless:
        """Returns the forward-scattering efficiency."""
        return self._cpp_Qforward * ureg.AU

    @property
    def Qratio(self) -> Dimensionless:
        """Returns the efficiency ratio of backscattering over total scattering."""
        return self._cpp_Qratio * ureg.AU

    @property
    def g(self) -> Dimensionless:
        """Returns the anisotropy factor."""
        return self._cpp_g * ureg.AU

    @property
    def Qpr(self) -> Dimensionless:
        """Returns the radiation pressure efficiency."""
        return self._cpp_Qpr * ureg.AU

    @property
    def Csca(self) -> Area:
        """Returns the scattering cross-section."""
        return (self._cpp_Csca * ureg.meter**2).to_compact()

    @property
    def Cext(self) -> Area:
        """Returns the extinction cross-section."""
        return (self._cpp_Cext * ureg.meter**2).to_compact()

    @property
    def Cabs(self) -> Area:
        """Returns the absorption cross-section."""
        return (self._cpp_Cabs * ureg.meter**2).to_compact()

    @property
    def Cpr(self) -> Area:
        """Returns the radiation pressure cross-section."""
        return (self._cpp_Cpr * ureg.meter**2).to_compact()

    @property
    def Cback(self) -> Area:
        """Returns the backscattering cross-section."""
        return (self._cpp_Cback * ureg.meter**2).to_compact()

    @property
    def Cforward(self) -> Area:
        """Returns the forward-scattering cross-section."""
        return (self._cpp_Cforward * ureg.meter**2).to_compact()

    @property
    def Cratio(self) -> Area:
        """Returns the ratio of backscattering cross-section over total scattering."""
        return (self._cpp_Cratio * ureg.meter**2).to_compact()

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

    def get_s1s2(
        self, sampling: int = 200, distance: Length = 1 * ureg.meter
    ) -> representations.S1S2:
        r"""
        Compute the S1 and S2 scattering amplitude functions for a spherical scatterer.

        The S1 and S2 parameters represent the scattering amplitudes for perpendicular and parallel polarizations of light, respectively. These parameters are fundamental in Mie theory, which describes the scattering of electromagnetic waves by spherical particles.

        The formulas for \( S_1 \) and \( S_2 \) are:

        .. math::
            S_1 = \sum\limits_{n=1}^{n_{\text{max}}} \frac{2n+1}{n(n+1)} \left( a_n \pi_n + b_n \tau_n \right) \\
            S_2 = \sum\limits_{n=1}^{n_{\text{max}}} \frac{2n+1}{n(n+1)} \left( a_n \tau_n + b_n \pi_n \right)

        Where:

        - :math:`a_n` and :math:`b_n`: Mie coefficients, which depend on the size, shape, and refractive index of the scatterer.
        - :math:`\pi_n` and :math:`\tau_n` \)`: Angular functions related to the angular components of the incident and scattered fields.
        - :math:`n_{\text{max}}`: Maximum number of terms in the series, determined by the size parameter of the scatterer.

        These scattering amplitude functions are essential for calculating properties such as scattering phase functions, efficiencies, and angular distribution of scattered light.

        Parameters
        ----------
        sampling : int
            The number of angular points used to sample the S1 and S2 functions. Higher sampling improves the resolution of the scattering pattern but increases computation time.
        distance : Length, optional
            The distance from the scatterer at which the S1 and S2 parameters are evaluated. This is typically set to 1 meter by default, but can be adjusted for specific setups.

        Returns
        -------
        representations.S1S2
            An object containing the computed S1 and S2 parameters, representing the scattering amplitudes for the two polarization components.

        Notes
        -----
        - The S1 and S2 parameters are central to Mie scattering theory and are used to derive many important scattering properties, such as intensity distributions and polarization effects.
        - The `sampling` parameter controls how finely the angular distribution is resolved. A higher value of `sampling` provides more detailed scattering information, which can be critical for accurately modeling the far-field pattern.

        Example
        -------
        You can use this method to compute the scattering properties of spherical particles, particularly in experiments where the polarization and scattering pattern of the light are important.

        Example usage:

        >>> s1s2 = scatterer.get_s1s2(sampling=500)
        >>> print(s1s2.S1, s1s2.S2)

        """
        return representations.S1S2(
            scatterer=self, sampling=sampling, distance=distance
        )

    def get_stokes(
        self, sampling: int = 200, distance: Length = 1 * ureg.meter
    ) -> representations.Stokes:
        r"""
        Compute the four Stokes parameters: I, Q, U, and V, which describe the polarization state of scattered light.

        The Stokes parameters provide a complete description of the polarization state of electromagnetic radiation, and are defined as follows:

        .. math::
            I &= \big| E_x \big|^2 + \big| E_y \big|^2 \quad &\text{(total intensity)} \\
            Q &= \big| E_x \big|^2 - \big| E_y \big|^2 \quad &\text{(linear polarization along x and y)} \\
            U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\} \quad &\text{(linear polarization at +45° and -45°)} \\
            V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\} \quad &\text{(circular polarization)}

        Where:

        - :math:`I`: Total intensity of the scattered light.
        - :math:`Q`: The degree of linear polarization along the x and y axes.
        - :math:`U`: The degree of linear polarization at +45° and -45° to the axes.
        - :math:`V`: The degree of circular polarization.

        These parameters provide a powerful way to represent the full polarization state of the scattered light.

        Parameters
        ----------
        sampling : int
            The number of angular points used to sample the Stokes parameters. A higher sampling value increases the resolution of the computed polarization pattern.
        distance : Length, optional
            The distance from the scatterer at which the Stokes parameters are evaluated. The default is 1 meter, but this can be adjusted based on the experimental setup.

        Returns
        -------
        Stokes
            An object containing the computed Stokes parameters (I, Q, U, V), which describe the polarization state of the scattered light.

        Notes
        -----
        - The Stokes parameters are essential for understanding and characterizing the polarization of scattered light. They are used in a wide variety of applications, including atmospheric optics, remote sensing, and optical communication.
        - The `sampling` parameter controls the angular resolution of the calculated Stokes parameters. Higher values provide finer detail, which can be important in accurately modeling polarization effects in scattering experiments.

        Example
        -------
        You can use this method to compute the Stokes parameters for a spherical scatterer, helping you analyze how the scatterer affects the polarization state of light.

        Example usage:

        >>> stokes_params = scatterer.get_stokes(sampling=500)
        >>> print(stokes_params.I, stokes_params.Q, stokes_params.U, stokes_params.V)

        """
        return representations.Stokes(
            scatterer=self, sampling=sampling, distance=distance
        )

    def get_farfield(
        self, sampling: int = 200, distance: Length = 1 * ureg.meter
    ) -> representations.FarField:
        r"""
        Compute the far-field scattering pattern for the scatterer.

        The far fields describe the behavior of the scattered electromagnetic waves at a large distance from the scatterer, where the waves can be approximated as planar.
        The computed far fields represent the scattered electric field components in the directions parallel and perpendicular to the plane of incidence.

        The far fields are computed as:

        .. math::
            \text{Fields} = E_{||}(\phi, \theta)^2, \, E_{\perp}(\phi, \theta)^2

        These components represent the intensities of the scattered electric field in the parallel (\( E_{||} \)) and perpendicular (\( E_{\perp} \)) directions as functions of the spherical angles \( \phi \) (azimuthal) and \( \theta \) (polar).

        The fields are expressed up to a constant phase factor:

        .. math::
            \exp{\left(-i k r \right)}

        Where:

        - :math:`k`: The wave number, related to the wavelength of the incident light.
        - :math:`r`: The distance from the scatterer (assumed to be large in the far-field approximation).

        Parameters
        ----------
        sampling : int
            The number of angular points used to sample the far-field pattern. A higher sampling value increases the angular resolution of the computed far-field intensities.
        distance : Length, optional
            The distance from the scatterer at which the far fields are evaluated. By default, this is set to 1 meter, but it can be adjusted to suit the specific experimental configuration.

        Returns
        -------
        representations.FarField
            An object containing the computed far-field components (parallel and perpendicular intensities), which represent the angular distribution of the scattered light in the far field.

        Notes
        -----
        - The far-field approximation assumes that the distance from the scatterer is large enough that the curvature of the wavefronts can be neglected, and the fields can be treated as planar.
        - This method is useful for determining the angular distribution of the scattered light, which is critical in applications such as radar cross-section analysis, optical scattering experiments, and remote sensing.

        Example
        -------
        You can use this method to compute the far-field scattering pattern of a spherical scatterer, which is essential in understanding how the scatterer redirects light in the far field.

        Example usage:

        >>> farfield = scatterer.get_farfield(sampling=500)
        >>> print(farfield.E_phi, farfield.E_theta)

        """
        return representations.FarField(
            scatterer=self, sampling=sampling, distance=distance
        )

    def get_spf(
        self, sampling: int = 200, distance: Length = 1 * ureg.meter
    ) -> representations.SPF:
        r"""
        Compute the scattering phase function (SPF) for the scatterer.

        The scattering phase function describes how the intensity of scattered light varies as a function of the scattering angles \( \theta \) (polar) and \( \phi \) (azimuthal). It is a key quantity in light scattering, as it characterizes the angular distribution of scattered light.

        The scattering phase function is computed as:

        .. math::
            \text{SPF} = \sqrt{ E_{\parallel}(\phi, \theta)^2 + E_{\perp}(\phi, \theta)^2 }

        Where:

        - :math:`E_{\parallel}`: The parallel component of the scattered electric field.
        - :math:`E_{\perp}`: The perpendicular component of the scattered electric field.
        - :math:`\phi` and :math:`\theta`: The azimuthal and polar angles, respectively, which describe the angular position of the scattered light.

        The SPF combines the intensity contributions from both polarization components to describe the total scattered intensity in different directions.

        Parameters
        ----------
        sampling : int
            The number of angular points used to sample the scattering phase function. A higher sampling value increases the angular resolution of the computed SPF.
        distance : Length, optional
            The distance from the scatterer at which the scattering phase function is evaluated. By default, this is set to 1 meter, but it can be adjusted depending on the specific experimental configuration.

        Returns
        -------
        representations.SPF
            An object containing the computed scattering phase function (SPF), representing the angular distribution of scattered light intensity.

        Notes
        -----
        - The scattering phase function is a critical quantity in the study of light scattering, as it provides insight into the directionality of the scattered light. It is commonly used in fields like atmospheric optics, biomedical imaging, and remote sensing.
        - The `sampling` parameter controls the angular resolution of the SPF. Higher sampling values yield more detailed scattering patterns, especially important when capturing small angular features.

        Example
        -------
        You can use this method to compute the scattering phase function of a spherical scatterer, helping to understand the angular distribution of light scattered by the particle.

        Example usage:

        >>> spf = scatterer.get_spf(sampling=500)
        >>> print(spf)

        """
        return representations.SPF(scatterer=self, sampling=sampling, distance=distance)

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

    def compute_nearfield(
        self, x: Length, y: Length, z: Length, radius: Length, field_type: str = "|E|"
    ):
        """
        Python wrapper for C++ near-field computation.

        This method bridges the gap between the C++ backend (_cpp_compute_near_field)
        and the Python NearField representation class.

        Parameters
        ----------
        x : Length
            Array of x coordinates of observation points.
        y : Length
            Array of y coordinates of observation points.
        z : Length
            Array of z coordinates of observation points.
        radius : Length
            Scatterer radius for inside/outside determination.
        field_type : str
            Field component type: "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "|E|", "|H|"

        Returns
        -------
        numpy.ndarray
            Complex array of field values at specified points.
        """
        return (
            self._cpp_compute_nearfields(
                x=x.to("meter").magnitude,
                y=y.to("meter").magnitude,
                z=z.to("meter").magnitude,
                field_type=field_type,
                radius=radius.to("meter").magnitude,
            )
            * ureg.volt
            / ureg.meter
        )

    def get_footprint(self, detector) -> representations.Footprint:
        r"""
        Compute the footprint of the scattered light coupling with the detector.

        The footprint represents the spatial distribution of the scattered light as it couples with the detector.
        This is important for understanding how the scattered light interacts with the detector's field of view and is influenced by both the scatterer's properties and the detector configuration.

        The footprint is computed using an inverse Fourier transform:

        .. math::
            \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi } (\xi, \nu),\
                \tilde{ \phi}_{l,m}(\xi, \nu)  \big\}
                (\delta_x, \delta_y) \big|^2

        Where:

        - :math:`\tilde{ \psi } (\xi, \nu)`: The Fourier transform of the scattered field.
        - :math:`\tilde{ \phi}_{l,m} (\xi, \nu)`: The Fourier transform of the detector's capturing field.
        - :math:`\mathscr{F}^{-1}`: The inverse Fourier transform operator.
        - :math:`(\delta_x, \delta_y)`: Spatial coordinates describing the footprint in the detector plane.

        The inverse Fourier transform yields the spatial pattern of the scattered light's interaction with the detector, which is then squared to compute the intensity of the footprint.

        Parameters
        ----------
        detector : GenericDetector
            The detector object that defines the capturing field and geometry. This object provides information about the detector's configuration, including its numerical aperture, position, and sensitivity.

        Returns
        -------
        representations.Footprint
            An object containing the computed scatterer footprint, representing the spatial distribution of the scattered light on the detector plane.

        Notes
        -----
        - The footprint provides valuable information about the coupling efficiency between the scattered light and the detector. This is useful in designing experiments where precise control of light detection is important, such as in microscopy, remote sensing, or optical measurements.
        - The interaction between the scattered field and the detector is influenced by factors such as the size, shape, and refractive index of the scatterer, as well as the detector's aperture and positioning.

        Example
        -------
        You can use this method to compute the footprint of the scattered light on a detector, helping to visualize how light from the scatterer is spatially distributed on the detector.

        Example usage:

        >>> footprint = scatterer.get_footprint(detector=my_detector)
        >>> print(footprint)

        """
        return representations.Footprint(scatterer=self, detector=detector)

    def _assign_index_or_material(
        self, property: RefractiveIndex | BaseMaterial
    ) -> tuple[RefractiveIndex | None, BaseMaterial | None]:
        """
        Determines whether the provided property is a refractive index (Quantity) or a material (BaseMaterial),
        and returns the corresponding values.

        Parameters:
        ----------
        property : RefractiveIndex or BaseMaterial
            The core property to be assigned, which can either be a refractive index (Quantity) or a material (BaseMaterial).

        Returns:
        -------
        tuple[RefractiveIndex | None, BaseMaterial | None]
            A tuple where the first element is the refractive index (Quantity) if provided, otherwise None.
            The second element is the material (BaseMaterial) if provided, otherwise None.

        Raises:
        ------
        ValueError:
            If the provided property is neither a Quantity (refractive index) nor a BaseMaterial.
        """
        if isinstance(property, RefractiveIndex):
            return property, None
        if isinstance(property, BaseMaterial):
            return (
                numpy.atleast_1d(
                    property.compute_refractive_index(
                        self.source.wavelength.to("meter").magnitude
                    )
                )[0]
                * ureg.RIU,
                property,
            )

        raise ValueError(
            f"Invalid material property: {property}. Expected a BaseMaterial or Quantity (RIU)."
        )
