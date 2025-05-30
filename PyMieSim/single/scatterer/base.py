#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Tuple
import numpy
from PyOptik.material.base_class import BaseMaterial
from tabulate import tabulate
from PyMieSim.single.representations import S1S2, FarField, Stokes, SPF, Footprint
from PyMieSim import units

class BaseScatterer(units.UnitsValidation):
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
    def size_parameter(self) -> units.Quantity:
        """Returns the size parameter of the scatterer."""
        return self._cpp_size_parameter * units.AU

    @property
    def cross_section(self) -> units.Quantity:
        """Returns the cross-section of the scatterer."""
        return (self._cpp_cross_section * units.meter**2).to_compact()

    @property
    def Qsca(self) -> units.Quantity:
        """Returns the scattering efficiency."""
        return self._cpp_Qsca * units.AU

    @property
    def Qext(self) -> units.Quantity:
        """Returns the extinction efficiency."""
        return self._cpp_Qext * units.AU

    @property
    def Qabs(self) -> units.Quantity:
        """Returns the absorption efficiency."""
        return self._cpp_Qabs * units.AU

    @property
    def Qback(self) -> units.Quantity:
        """Returns the backscattering efficiency."""
        return self._cpp_Qback * units.AU

    @property
    def Qforward(self) -> units.Quantity:
        """Returns the forward-scattering efficiency."""
        return self._cpp_Qforward * units.AU

    @property
    def Qratio(self) -> units.Quantity:
        """Returns the efficiency ratio of backscattering over total scattering."""
        return self._cpp_Qratio * units.AU

    @property
    def g(self) -> units.Quantity:
        """Returns the anisotropy factor."""
        return self._cpp_g * units.AU

    @property
    def Qpr(self) -> units.Quantity:
        """Returns the radiation pressure efficiency."""
        return self._cpp_Qpr * units.AU

    @property
    def Csca(self) -> units.Quantity:
        """Returns the scattering cross-section."""
        return (self._cpp_Csca * units.meter ** 2).to_compact()

    @property
    def Cext(self) -> units.Quantity:
        """Returns the extinction cross-section."""
        return (self._cpp_Cext * units.meter ** 2).to_compact()

    @property
    def Cabs(self) -> units.Quantity:
        """Returns the absorption cross-section."""
        return (self._cpp_Cabs * units.meter ** 2).to_compact()

    @property
    def Cpr(self) -> units.Quantity:
        """Returns the radiation pressure cross-section."""
        return (self._cpp_Cpr * units.meter ** 2).to_compact()

    @property
    def Cback(self) -> units.Quantity:
        """Returns the backscattering cross-section."""
        return (self._cpp_Cback * units.meter ** 2).to_compact()

    @property
    def Cforward(self) -> units.Quantity:
        """Returns the forward-scattering cross-section."""
        return (self._cpp_Cforward * units.meter ** 2).to_compact()

    @property
    def Cratio(self) -> units.Quantity:
        """Returns the ratio of backscattering cross-section over total scattering."""
        return (self._cpp_Cratio * units.meter ** 2).to_compact()

    def get_farfields_array(self, phi: numpy.ndarray, theta: numpy.ndarray, r: units.Quantity) -> Tuple[numpy.ndarray, numpy.ndarray]:
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
        return self._cpp_get_fields(phi=phi, theta=theta, r=r.to_base_units().magnitude)

    def get_s1s2(self, sampling: int = 200, distance: units.Quantity = 1 * units.meter) -> S1S2:
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
        distance : units.Quantity, optional
            The distance from the scatterer at which the S1 and S2 parameters are evaluated. This is typically set to 1 meter by default, but can be adjusted for specific setups.

        Returns
        -------
        S1S2
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
        return S1S2(scatterer=self, sampling=sampling, distance=distance)

    def get_stokes(self, sampling: int = 200, distance: units.Quantity = 1 * units.meter) -> Stokes:
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
        distance : units.Quantity, optional
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
        return Stokes(scatterer=self, sampling=sampling, distance=distance)

    def get_far_field(self, sampling: int = 200, distance: units.Quantity = 1 * units.meter) -> FarField:
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
        distance : units.Quantity, optional
            The distance from the scatterer at which the far fields are evaluated. By default, this is set to 1 meter, but it can be adjusted to suit the specific experimental configuration.

        Returns
        -------
        FarField
            An object containing the computed far-field components (parallel and perpendicular intensities), which represent the angular distribution of the scattered light in the far field.

        Notes
        -----
        - The far-field approximation assumes that the distance from the scatterer is large enough that the curvature of the wavefronts can be neglected, and the fields can be treated as planar.
        - This method is useful for determining the angular distribution of the scattered light, which is critical in applications such as radar cross-section analysis, optical scattering experiments, and remote sensing.

        Example
        -------
        You can use this method to compute the far-field scattering pattern of a spherical scatterer, which is essential in understanding how the scatterer redirects light in the far field.

        Example usage:

        >>> far_field = scatterer.get_far_field(sampling=500)
        >>> print(far_field.E_phi, far_field.E_theta)

        """
        return FarField(scatterer=self, sampling=sampling, distance=distance)

    def get_spf(self, sampling: int = 200, distance: units.Quantity = 1 * units.meter) -> SPF:
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
        distance : units.Quantity, optional
            The distance from the scatterer at which the scattering phase function is evaluated. By default, this is set to 1 meter, but it can be adjusted depending on the specific experimental configuration.

        Returns
        -------
        SPF
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
        return SPF(scatterer=self, sampling=sampling, distance=distance)

    def get_footprint(self, detector) -> Footprint:
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
        Footprint
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
        return Footprint(scatterer=self, detector=detector)

    def _assign_index_or_material(self, property: units.Quantity | BaseMaterial) -> tuple[units.Quantity | None, BaseMaterial | None]:
        """
        Determines whether the provided property is a refractive index (Quantity) or a material (BaseMaterial),
        and returns the corresponding values.

        Parameters:
        ----------
        property : units.Quantity or BaseMaterial
            The core property to be assigned, which can either be a refractive index (Quantity) or a material (BaseMaterial).

        Returns:
        -------
        tuple[units.Quantity | None, BaseMaterial | None]
            A tuple where the first element is the refractive index (Quantity) if provided, otherwise None.
            The second element is the material (BaseMaterial) if provided, otherwise None.

        Raises:
        ------
        ValueError:
            If the provided property is neither a Quantity (refractive index) nor a BaseMaterial.
        """
        if isinstance(property, units.Quantity):
            return property, None
        if isinstance(property, BaseMaterial):
            return numpy.atleast_1d(property.compute_refractive_index(self.source.wavelength.to_base_units().magnitude))[0] * units.RIU, property

        raise ValueError(f"Invalid material property: {property}. Expected a BaseMaterial or Quantity (RIU).")
