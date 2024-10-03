#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from PyOptik.base_class import BaseMaterial
from pydantic.dataclasses import dataclass
from tabulate import tabulate
from PyMieSim.single.source.base import BaseSource
from PyMieSim.single.representations import S1S2, FarField, Stokes, SPF, Footprint
from PyMieSim.units import RIU, Quantity, meter, AU

config_dict = dict(
    arbitrary_types_allowed=True,
    kw_only=True,
    slots=True,
    extra='forbid'
)

@dataclass(config=config_dict, kw_only=True)
class BaseScatterer:
    """
    A generic class for a scatterer, providing properties and methods to compute various scattering-related quantities.

    Attributes:
        source (Union[source.PlaneWave, source.Gaussian]): Light source object containing info on polarization and wavelength.
        medium_index (float): Refractive index of scatterer medium. Default is 1.0.
        material (Optional[BaseMaterial]): BaseMaterial of which the scatterer is made, if index is not specified. Default is None.
    """
    source: BaseSource
    medium_property: Quantity | BaseMaterial

    def __post_init__(self):
        self.medium_index, self.medium_material = self._assign_index_or_material(self.medium_property)

        self.set_binding()

    def print_properties(self) -> None:
        """Prints a table of the scatterer's properties and their values."""
        data = [getattr(self, name) for name in self.property_names]
        property_dict = {"Property": self.property_names, "Value": data}

        table = tabulate(property_dict, headers="keys")
        print(table)

    @property
    def size_parameter(self) -> Quantity:
        """Returns the size parameter of the scatterer."""
        return self.binding.size_parameter * AU

    @property
    def area(self) -> Quantity:
        """Returns the area of the scatterer."""
        return (self.binding.area * meter**2).to_compact()

    @property
    def Qsca(self) -> Quantity:
        """Returns the scattering efficiency."""
        return self.binding.Qsca * AU

    @property
    def Qext(self) -> Quantity:
        """Returns the extinction efficiency."""
        return self.binding.Qext * AU

    @property
    def Qabs(self) -> Quantity:
        """Returns the absorption efficiency."""
        return self.binding.Qabs * AU

    @property
    def Qback(self) -> Quantity:
        """Returns the backscattering efficiency."""
        return self.binding.Qback * AU

    @property
    def Qforward(self) -> Quantity:
        """Returns the forward-scattering efficiency."""
        return self.binding.Qforward * AU

    @property
    def Qratio(self) -> Quantity:
        """Returns the efficiency ratio of backscattering over total scattering."""
        return self.binding.Qratio * AU

    @property
    def g(self) -> Quantity:
        """Returns the anisotropy factor."""
        return self.binding.g * AU

    @property
    def Qpr(self) -> Quantity:
        """Returns the radiation pressure efficiency."""
        return self.binding.Qpr * AU

    @property
    def Csca(self) -> Quantity:
        """Returns the scattering cross-section."""
        return (self.binding.Csca * meter ** 2).to_compact()

    @property
    def Cext(self) -> Quantity:
        """Returns the extinction cross-section."""
        return (self.binding.Cext * meter ** 2).to_compact()

    @property
    def Cabs(self) -> Quantity:
        """Returns the absorption cross-section."""
        return (self.binding.Cabs * meter ** 2).to_compact()

    @property
    def Cpr(self) -> Quantity:
        """Returns the radiation pressure cross-section."""
        return (self.binding.Cpr * meter ** 2).to_compact()

    @property
    def Cback(self) -> Quantity:
        """Returns the backscattering cross-section."""
        return (self.binding.Cback * meter ** 2).to_compact()

    @property
    def Cforward(self) -> Quantity:
        """Returns the forward-scattering cross-section."""
        return (self.binding.Cforward * meter ** 2).to_compact()

    @property
    def Cratio(self) -> Quantity:
        """Returns the ratio of backscattering cross-section over total scattering."""
        return (self.binding.Cratio * meter ** 2).to_compact()

    def get_farfields_array(self, phi: numpy.ndarray, theta: numpy.ndarray, r: numpy.ndarray) -> numpy.ndarray:
        """
        Computes the scattering far field for unstructured coordinates.

        The method computes the fields up to a constant phase value.

        Args:
            phi (numpy.ndarray): The phi angles in radians.
            theta (numpy.ndarray): The theta angles in radians.
            r (numpy.ndarray): The radial distances.

        Returns:
            numpy.ndarray: The computed far fields.
        """
        return self.binding.get_fields(phi=phi, theta=theta, r=r)

    def get_s1s2(self, **kwargs) -> S1S2:
        r"""
        Computes the S1 and S2 parameters.

        These parameters are computed for spherical scatterers using the following formulas:

        .. math::
            S_1=\sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)}(a_n \pi_n+b_n \tau_n) \\
            S_2=\sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)}(a_n \tau_n+b_n \pi_n)

        Args:
            kwargs Additional keyword arguments.

        Returns:
            S1S2: The computed S1 and S2 parameters.
        """
        return S1S2(scatterer=self, **kwargs)

    def get_stokes(self, **kwargs) -> Stokes:
        r"""
        Computes the four Stokes components: I, Q, U, and V.

        These components are defined as:

        .. math::
            I &= \big| E_x \big|^2 + \big| E_y \big|^2 \\
            Q &= \big| E_x \big|^2 - \big| E_y \big|^2 \\
            U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\} \\
            V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\}

        Args:
            kwargs Additional keyword arguments.

        Returns:
            Stokes: The computed Stokes parameters.
        """
        return Stokes(scatterer=self, **kwargs)

    def get_far_field(self, **kwargs) -> FarField:
        r"""
        Computes the scattering far fields.

        The far fields are computed as:

        .. math::
            \text{Fields} = E_{||}(\phi,\theta)^2, E_{\perp}(\phi,\theta)^2

        The fields are up to a constant phase value:

        .. math::
            \exp{\big(-i k r \big)}

        Args:
            kwargs Additional keyword arguments.

        Returns:
            FarField: The computed far fields.
        """
        return FarField(scatterer=self, **kwargs)

    def get_spf(self, **kwargs) -> SPF:
        r"""
        Computes the scattering phase function.

        The scattering phase function is computed as:

        .. math::
            \text{SPF} = \sqrt{ E_{\parallel}(\phi,\theta)^2 + E_{\perp}(\phi,\theta)^2 }

        Args:
            kwargs Additional keyword arguments.

        Returns:
            SPF: The computed scattering phase function.
        """
        return SPF(scatterer=self, **kwargs)

    def get_footprint(self, detector) -> Footprint:
        r"""
        Computes the footprint of the scattered light coupling with the detector.

        The footprint is computed as:

        .. math::
            \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi } (\xi, \nu),\
                   \tilde{ \phi}_{l,m}(\xi, \nu)  \big\}
                   (\delta_x, \delta_y) \big|^2

        Args:
            detector (GenericDetector): The detector object.

        Returns:
            Footprint: The computed scatterer footprint.
        """
        return Footprint(scatterer=self, detector=detector)

    def _assign_index_or_material(self, property: Quantity | BaseMaterial) -> tuple[Quantity | None, BaseMaterial | None]:
        """
        Determines whether the provided property is a refractive index (Quantity) or a material (BaseMaterial),
        and returns the corresponding values.

        Parameters:
        ----------
        property : Quantity or BaseMaterial
            The core property to be assigned, which can either be a refractive index (Quantity) or a material (BaseMaterial).

        Returns:
        -------
        tuple[Quantity | None, BaseMaterial | None]
            A tuple where the first element is the refractive index (Quantity) if provided, otherwise None.
            The second element is the material (BaseMaterial) if provided, otherwise None.

        Raises:
        ------
        ValueError:
            If the provided property is neither a Quantity (refractive index) nor a BaseMaterial.
        """
        if isinstance(property, Quantity):
            return property, None
        if isinstance(property, BaseMaterial):
            return numpy.atleast_1d(property.compute_refractive_index(self.source.wavelength.to_base_units().magnitude))[0] * RIU, property
        raise ValueError(f"Invalid material property: {property}. Expected a BaseMaterial or Quantity (RIU).")
