#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import MPSPlots
from pydantic.dataclasses import dataclass
from PyMieSim.special_functions import spherical_to_cartesian, rotate_on_x
from typing import List
import pyvista
from MPSPlots.colormaps import blue_black_red
import matplotlib.pyplot as plt
from PyMieSim.units import Quantity, meter


config_dict = dict(
    arbitrary_types_allowed=True,
    kw_only=True,
    slots=True,
    extra='forbid'
)


@dataclass(config=config_dict, kw_only=True)
class BaseRepresentation():
    """
    Base class for scattering representations.

    Parameters
    ----------
    scatterer : BaseScatterer
        The scatterer object, representing the physical scatterer in the simulation.
    sampling : int
        The number of points used for evaluating the Stokes parameters in spherical coordinates (default is 100).
    distance : float
        The distance from the scatterer at which fields are evaluated (default is 1.0).

    Methods:
        compute_components: A placeholder method intended to be overridden by subclasses for computing specific scattering components.
    """
    scatterer: object
    sampling: int
    distance: Quantity

    def __post_init__(self):
        fields = self.scatterer._cpp_get_full_fields(
            sampling=self.sampling,
            distance=self.distance.to_base_units().magnitude
        )

        self.E_phi, self.E_theta, self.theta, self.phi = fields

        self.compute_components()

    def compute_components(self) -> None:
        """
        Placeholder method for computing scattering components. Intended to be overridden by subclasses.
        """
        raise NotImplementedError("This method should be implemented by subclasses.")

    def get_colormap_limits(self, scalar: numpy.ndarray, symmetric: bool = False):
        if symmetric:
            max_abs = numpy.abs(scalar).max()
            return [-max_abs, max_abs]
        else:
            return None

    def add_theta_vector_to_3d_plot(
            self,
            scene: pyvista.Plotter,
            n_points: int = 20,
            opacity: float = 1.0,
            radius: float = 1.0,
            color: str = 'black') -> None:
        """
        Adds a vector field to the 3D plot, representing vectors in the theta direction.

        Parameters
        ----------
        scene : pyvista.Plotter
            The 3D plotting scene to which the vectors will be added.
        n_points : int
            Number of points to generate along the theta and phi directions. Default is 100.
        opacity : float
            Opacity of the vectors. Default is 1.0.
        radius : float
            Radius at which to place the vectors. Default is 1.0.
        color : str
            Color of the vectors. Default is 'black'.
        """
        theta = numpy.linspace(0, 360, n_points)
        phi = numpy.linspace(180, 0, n_points)

        # Define the vector direction (unit vector along x-axis)
        vector = numpy.array([1, 0, 0])

        # Convert spherical coordinates to Cartesian coordinates
        x, y, z = pyvista.transform_vectors_sph_to_cart(theta, phi, radius, *vector)

        # Combine the Cartesian coordinates into a vector array
        vector_field = numpy.c_[x.ravel(), y.ravel(), z.ravel()]

        # Create a structured grid from spherical coordinates
        spherical_grid = pyvista.grid_from_sph_coords(theta, phi, radius)
        spherical_grid.point_data["component"] = vector_field * 0.1

        # Generate glyphs (arrows) for the vectors
        glyphs = spherical_grid.glyph(orient="component", scale="component", tolerance=0.005)

        # Add the vector glyphs to the scene
        scene.add_mesh(glyphs, color=color, opacity=opacity)

    def add_phi_vector_to_3d_plot(
        self,
        scene: pyvista.Plotter,
        n_points: int = 20,
        opacity: float = 1.0,
        radius: float = 1.0,
        color: str = 'black') -> None:
        """
        Adds a vector field to the 3D plot, representing vectors in the phi direction.

        Parameters
        ----------
        scene : pyvista.Plotter
            The 3D plotting scene to which the vectors will be added.
        n_points : int
            Number of points to generate along the theta and phi directions. Default is 100.
        opacity : float
            Opacity of the vectors. Default is 1.0.
        radius : float
            Radius at which to place the vectors. Default is 1.0.
        color : str
            Color of the vectors. Default is 'black'.
        """
        theta = numpy.linspace(0, 360, n_points)
        phi = numpy.linspace(180, 0, n_points)

        # Define the vector direction (unit vector along y-axis)
        vector = numpy.array([0, 1, 0])

        # Convert spherical coordinates to Cartesian coordinates
        x, y, z = pyvista.transform_vectors_sph_to_cart(theta, phi, radius, *vector)

        # Combine the Cartesian coordinates into a vector array
        vector_field = numpy.c_[x.ravel(), y.ravel(), z.ravel()]

        # Create a structured grid from spherical coordinates
        spherical_grid = pyvista.grid_from_sph_coords(theta, phi, radius)
        spherical_grid.point_data["component"] = vector_field * 0.1

        # Generate glyphs (arrows) for the vectors
        glyphs = spherical_grid.glyph(orient="component", scale="component", tolerance=0.005)

        # Add the vector glyphs to the scene
        scene.add_mesh(glyphs, color=color, opacity=opacity)


@dataclass(config=config_dict, kw_only=True)
class Stokes(BaseRepresentation):
    r"""
    Represents the scattering far-field in the Stokes representation.

    Inherits from BaseRepresentation and calculates the Stokes parameters which describe the polarization state of light.

    The Stokes parameters (I, Q, U, V) are defined according to their conventional definitions, representing the total intensity,
    difference in intensities between horizontal and vertical polarizations, difference in intensities between two diagonal polarizations,
    and the right and left circular polarizations, respectively.

    | The stokes parameters are:
    |     I : intensity of the fields
    |     Q : linear polarization parallel to incident polarization
    |     U : linear polarization 45 degree to incident polarization
    |     V : Circular polarization

    .. math:
        I &= \big| E_x \big|^2 + \big| E_y \big|^2 \\[10pt]

        Q &= \big| E_x \big|^2 - \big| E_y \big|^2 \\[10pt]

        U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\} \\[10pt]

        V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\} \\[10pt]

    Methods:
        compute_components: Computes the Stokes parameters based on the electric field components.
        plot: Visualizes the Stokes parameters on a 3D plot.

    """

    def compute_components(self) -> None:
        r"""
        Computes the Stokes parameters (I, Q, U, V) based on the electric field components (E_phi and E_theta).

        The method calculates the normalized intensity (I), linear polarizations (Q and U), and circular polarization (V) of the light
        scattered by the particle, using the electric field components in spherical coordinates.

        The Stokes parameters are calculated using the following formulas:

        .. math:
            - I = |E_phi|^2 + |E_theta|^2
            - Q = |E_phi|^2 - |E_theta|^2
            - U = 2 * Re{E_phi * E_theta*}
            - V = -2 * Im{E_phi * E_theta*}

        The results are stored as attributes of the instance: I, Q, U, and V.

        """
        intensity = numpy.abs(self.E_phi)**2 + numpy.abs(self.E_theta)**2

        self.I = intensity / numpy.max(intensity)  # noqa: E741
        self.Q = (numpy.abs(self.E_phi)**2 - numpy.abs(self.E_theta)**2) / intensity
        self.U = (+2 * numpy.real(self.E_phi * self.E_theta.conjugate())) / intensity
        self.V = (-2 * numpy.imag(self.E_phi * self.E_theta.conjugate())) / intensity

    def plot(
        self,
        unit_size: List[float] = (400, 400),
        background_color: str = 'white',
        show_edges: bool = False,
        colormap: str = blue_black_red,
        opacity: float = 1.0,
        symmetric_colormap: bool = False,
        show_axis_label: bool = False) -> None:
        """
        Visualizes the Stokes parameters (I, Q, U, V) on a 3D plot.

        Parameters
        ----------
        unit_size : List[float]
            The size of each subplot in pixels (width, height). Default is (400, 400).
        background_color : str
            The background color of the plot. Default is 'white'.
        show_edges : bool
            If True, displays the edges of the mesh. Default is False.
        colormap : str
            The colormap to use for scalar mapping. Default is 'blue_black_red'.
        opacity : float
            The opacity of the mesh. Default is 1.0.
        symmetric_colormap : bool
            If True, the colormap will be symmetric around zero. Default is False.
        show_axis_label : bool
            If True, shows the axis labels. Default is False.
        """
        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)
        x, y, z = spherical_to_cartesian(r=numpy.full_like(phi_mesh, 0.5), phi=phi_mesh, theta=theta_mesh)

        window_size = (unit_size[1] * 4, unit_size[0])  # Four subplots horizontally

        scene = pyvista.Plotter(theme=pyvista.themes.DocumentTheme(), window_size=window_size, shape=(1, 4))
        scene.set_background(background_color)

        for idx, (name, field) in enumerate(zip(['I', 'Q', 'U', 'V'], [self.I, self.Q, self.U, self.V])):
            field = field.flatten(order='F')
            mesh = pyvista.StructuredGrid(x, y, z)
            scene.subplot(0, idx)

            colormap_limits = self.get_colormap_limits(
                scalar=field,
                symmetric=symmetric_colormap
            )

            mapping = scene.add_mesh(
                mesh,
                cmap=colormap,
                scalars=field, opacity=opacity,
                style='surface',
                show_edges=show_edges,
                clim=colormap_limits,
                show_scalar_bar=False
            )

            scene.add_axes_at_origin(labels_off=not show_axis_label)
            scene.add_scalar_bar(mapper=mapping.mapper, title=f'{name} field')

        scene.show()


@dataclass(config=config_dict, kw_only=True)
class FarField(BaseRepresentation):
    r"""
    Represents the scattering far-field in spherical coordinates.

    Inherits from BaseRepresentation and visualizes the far-field pattern characterized by the perpendicular and parallel components
    of the electric field in spherical coordinates.

    .. math::
        \text{Fields} = E_{||}(\phi,\theta)^2, E_{\perp}(\phi,\theta)^2

    Methods:
        compute_components: Calculates the field components. This implementation is a placeholder, as the components are precomputed.
        plot: Visualizes the far-field pattern in a 3D plot.

    """

    def compute_components(self) -> None:
        """
        Placeholder method in FarField class. Does not perform any computation as field components are precomputed.

        This method is intended to be consistent with the structure of BaseRepresentation but does not need to modify or compute
        any attributes for FarField instances.
        """
        return

    def plot(
            self,
            unit_size: List[float] = (400, 400),
            background_color: str = 'white',
            show_edges: bool = False,
            colormap: str = blue_black_red,
            opacity: float = 1.0,
            symmetric_colormap: bool = False,
            show_axis_label: bool = False) -> None:
        """
        Visualizes the Far field (in phi and theta vector projections) on a 3D plot.

        Parameters
        ----------
        unit_size : List[float]
            The size of each subplot in pixels (width, height). Default is (400, 400).
        background_color : str
            The background color of the plot. Default is 'white'.
        show_edges : bool
            If True, displays the edges of the mesh. Default is False.
        colormap : str
            The colormap to use for scalar mapping. Default is 'blue_black_red'.
        opacity : float
            The opacity of the mesh. Default is 1.0.
        symmetric_colormap : bool
            If True, the colormap will be symmetric around zero. Default is False.
        show_axis_label : bool
            If True, shows the axis labels. Default is False.
        """
        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)
        x, y, z = spherical_to_cartesian(r=numpy.full_like(phi_mesh, 0.5), phi=phi_mesh, theta=theta_mesh)

        window_size = (unit_size[1] * 4, unit_size[0])  # Two subplots horizontally

        scene = pyvista.Plotter(theme=pyvista.themes.DocumentTheme(), window_size=window_size, shape=(1, 4))
        scene.set_background(background_color)

        repr = [self.E_phi.real, self.E_phi.imag, self.E_theta.real, self.E_theta.imag]
        repr_label = ['phi real', 'phi imag', 'theta real', 'theta imag']

        for idx, (label, field) in enumerate(zip(repr_label, repr)):
            field = field.flatten(order='F')
            mesh = pyvista.StructuredGrid(x, y, z)
            scene.subplot(0, idx)

            colormap_limits = self.get_colormap_limits(
                scalar=field,
                symmetric=symmetric_colormap
            )

            mapping = scene.add_mesh(
                mesh,
                cmap=colormap,
                scalars=field,
                opacity=opacity,
                style='surface',
                show_edges=show_edges,
                clim=colormap_limits,
                show_scalar_bar=False
            )
            if 'theta' in label:
                self.add_theta_vector_to_3d_plot(scene=scene, radius=0.6)

            if 'phi' in label:
                self.add_phi_vector_to_3d_plot(scene=scene, radius=0.6)

            scene.add_axes_at_origin(labels_off=not show_axis_label)
            scene.add_scalar_bar(mapper=mapping.mapper, title=f'{label} field')

        scene.show()


@dataclass(config=config_dict, kw_only=True)
class SPF(BaseRepresentation):
    r"""
    Represents the Scattering Phase Function (SPF).

    Inherits from BaseRepresentation and computes the SPF, which is a measure of how light is scattered by a particle at different angles.

    .. math::
        \text{SPF} = E_{\parallel}(\phi,\theta)^2 + E_{\perp}(\phi,\theta)^2

    Methods:
        compute_components: Computes the SPF based on the electric field components.
        plot: Visualizes the SPF on a 3D plot.

    """

    def compute_components(self) -> None:
        """
        Computes the Scattering Phase Function (SPF) based on the electric field components (E_phi and E_theta).

        The SPF is calculated as the square root of the sum of the squared magnitudes of the electric field components, representing
        the total scattering intensity distribution as a function of angles.

        The result is stored as the SPF attribute of the instance.
        """
        self.SPF = numpy.sqrt(numpy.abs(self.E_phi)**2 + numpy.abs(self.E_theta)**2)

    def plot(
            self,
            unit_size: List[float] = (400, 400),
            background_color: str = 'white',
            show_edges: bool = False,
            colormap: str = 'viridis',
            opacity: float = 1.0,
            set_surface: bool = True,
            show_axis_label: bool = False) -> None:
        """
        Visualizes the scattering phase function on a 3D plot.

        This method creates a 3D visualization of the scattering phase function (SPF). It allows customization
        of the plot's appearance, including the colormap, mesh opacity, and whether or not to display mesh edges
        and axis labels.

        Parameters
        ----------
        unit_size : List[float]
            The size of the plot window in pixels (width, height). Default is (400, 400).
        background_color : str
            The background color of the plot. Default is 'white'.
        show_edges : bool
            If True, displays the edges of the mesh. Default is False.
        colormap : str
            The colormap to use for scalar mapping. Default is 'viridis'.
        opacity : float
            The opacity of the mesh. Default is 1.0.
        set_surface : bool
            If True, the surface represents the scaled SPF; if False, a unit sphere is used. Default is True.
        show_axis_label : bool
            If True, shows the axis labels. Default is False.
        """
        # Define the window size based on the unit size provided
        window_size = (unit_size[1], unit_size[0])  # One subplot

        # Create a PyVista plotting scene with the specified theme and window size
        scene = pyvista.Plotter(theme=pyvista.themes.DocumentTheme(), window_size=window_size)

        # Set the background color of the scene
        scene.set_background(background_color)

        # Add the 3D axis-aligned plot to the scene using the specified settings
        mapping = self._add_to_3d_ax(
            scene=scene,
            colormap=colormap,
            opacity=opacity,
            show_edges=show_edges,
            set_surface=set_surface
        )

        # Optionally add axis labels
        scene.add_axes_at_origin(labels_off=not show_axis_label)

        # Add a scalar bar to the scene to represent the scattering phase function
        scene.add_scalar_bar(mapper=mapping.mapper, title='Scattering Phase Function')

        # Display the scene
        scene.show()

    def _add_to_3d_ax(self, scene: pyvista.Plotter, set_surface: bool = False, show_edges: bool = False, colormap: str = 'viridis', opacity: float = 1.0) -> None:
        """
        Adds a 3D surface plot to the provided PyVista scene based on the scattering phase function (SPF).

        This method generates a 3D surface plot of the SPF using spherical coordinates, and adds it to the scene.
        The surface can either represent the actual SPF or a normalized unit sphere, depending on the `set_surface` flag.
        The appearance of the surface can be customized using various parameters.

        Parameters
        ----------
        scene : pyvista.Plotter
            The PyVista plotting scene where the surface will be added.
        set_surface : bool
            If True, the surface will represent the scaled SPF; if False, a unit sphere is used. Default is True.
        show_edges : bool
            If True, edges of the mesh will be displayed. Default is False.
        colormap : str
            The colormap to use for visualizing the scalar field. Default is 'viridis'.
        opacity : float
            The opacity of the surface mesh. Default is 1.0.
        """
        # Create mesh grids for phi and theta
        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)

        # Normalize the scattering phase function (SPF) for visualization
        scalar = self.SPF / self.SPF.max() * 2

        # Determine the coordinates based on whether the surface represents the SPF or a unit sphere
        if set_surface:
            x, y, z = spherical_to_cartesian(r=scalar, phi=phi_mesh, theta=theta_mesh)
        else:
            x, y, z = spherical_to_cartesian(r=numpy.ones(phi_mesh.shape) * 0.5, phi=phi_mesh, theta=theta_mesh)

        # Create a structured grid from the calculated coordinates
        mesh = pyvista.StructuredGrid(x, y, z)

        # Add the surface mesh to the scene
        mapping = scene.add_mesh(
            mesh,
            cmap=colormap,
            scalars=scalar.flatten(order='F'),
            opacity=opacity,
            style='surface',
            show_edges=show_edges,
            show_scalar_bar=False
        )

        return mapping


@dataclass(config=config_dict, kw_only=True)
class S1S2(BaseRepresentation):
    """
    Represents the S1 and S2 scattering functions, which are components of the scattering matrix.

    Parameters
    ----------
    scatterer : BaseScatterer
        The scatterer object.
    sampling : int
        Number of points for evaluating the S1 and S2 functions.
    distance : Quantity
        Distance at which the fields are evaluated.

    Methods:
        compute_components: Computes the S1 and S2 functions based on the scatterer's properties.
        plot: Visualizes the S1 and S2 functions on a polar plot.
    """
    def __post_init__(self):
        self.phi = numpy.linspace(-180, 180, self.sampling)

        self.compute_components()

    def compute_components(self) -> None:
        """
        Computes the S1 and S2 scattering parameters based on the scatterer's properties and the scattering angle phi.

        S1 and S2 are integral parts of the scattering matrix describing the change in polarization state of light upon scattering.

        The method calculates these parameters for a range of phi angles and stores them as the S1 and S2 attributes of the instance.
        """
        self.S1, self.S2 = self.scatterer._cpp_get_s1s2(phi=numpy.deg2rad(self.phi) + numpy.pi / 2)

    def plot(self) -> None:
        """
        Plots the S1 and S2 Stokes parameters on polar plots.

        The method generates two polar plots: one for the absolute values of the S1 parameter and another
        for the S2 parameter, filling the area between the radial axis and the parameter values.

        Returns
        -------
        None
            This method does not return a value. It displays the polar plots.
        """
        with plt.style.context(MPSPlots.styles.mps):
            figure, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'polar': True})

            # Plot for S1 parameter
            axes[0].set(title=r'S$_1$ parameter')
            axes[0].fill_between(
                numpy.deg2rad(self.phi),
                y1=0,
                y2=numpy.abs(self.S1),
                color='C0',
                edgecolor='black'
            )

            # Plot for S2 parameter
            axes[1].set(title=r'S$_2$ parameter')
            axes[1].fill_between(
                numpy.deg2rad(self.phi),
                y1=0,
                y2=numpy.abs(self.S2),
                color='C1',
                edgecolor='black'
            )

            plt.show()


@dataclass(config=config_dict, kw_only=True)
class Footprint():
    r"""
    Represents the footprint of the scatterer as detected by various detectors.

    .. math::
        \text{Footprint} = \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi }\
        (\xi, \nu), \tilde{ \phi}_{l,m}(\xi, \nu)  \big\} \
        (\delta_x, \delta_y) \big|^2

    Parameters
    ----------
    detector : BaseDetector
        The detector object.
    scatterer : BaseScatterer
        The scatterer object.
    sampling : int
        Number of points to evaluate the Stokes parameters in spherical coordinates (default is 500).
    padding_factor : int
        Padding factor for the Fourier transform (default is 20).

    Methods:
        compute_footprint: Computes the footprint based on the far-field patterns and detector characteristics.
        plot: Visualizes the computed footprint.
    """
    detector: object
    scatterer: object
    sampling: int = 200
    padding_factor: int = 20

    def __post_init__(self):
        self.compute_footprint()

    def compute_footprint(self):
        """
        Computes the footprint of the scatterer as detected by the specified detector.

        The footprint is calculated based on the far-field scattering patterns and the characteristics of the detector,
        using a Fourier transform to project the far-field onto the detector plane.

        The computed footprint and the corresponding spatial coordinates are stored as attributes of the instance.

        Warning: this function do not currently take account of the cache block on the detector.
        """
        max_angle = self.detector.max_angle
        n_point = complex(self.sampling)

        phi, theta = numpy.mgrid[
            -max_angle.to('radian').magnitude: max_angle.to('radian').magnitude: n_point, 0: numpy.pi: n_point
        ]

        max_distance_direct_space = 1 / (numpy.sin(max_angle) * self.scatterer.source.wavenumber / (2 * numpy.pi))

        x = y = numpy.linspace(-1, 1, self.sampling) * self.sampling / 2 * max_distance_direct_space / self.padding_factor

        _, phi, theta = rotate_on_x(phi + numpy.pi / 2, theta, numpy.pi / 2)

        far_field_para, far_field_perp = self.scatterer.get_farfields_array(
            phi=phi.ravel() + numpy.pi / 2,
            theta=theta.ravel(),
            r=1.0 * meter,
        )

        detector_structured_farfield = self.detector.get_structured_scalarfield(sampling=self.sampling)

        perpendicular_projection = detector_structured_farfield * far_field_perp.reshape(theta.shape)

        parallel_projection = detector_structured_farfield * far_field_para.reshape(theta.shape)

        fourier_parallel = self.get_fourier_component(parallel_projection)
        fourier_perpendicular = self.get_fourier_component(perpendicular_projection)

        self.mapping = (fourier_parallel + fourier_perpendicular)
        self.direct_x = x
        self.direct_y = y

    def get_fourier_component(self, scalar: numpy.ndarray) -> numpy.ndarray:
        """
        Computes the Fourier component of a given scalar field.

        This method performs a two-dimensional inverse Fourier transform on the input scalar field, which represents
        a projection (either parallel or perpendicular) of the far-field pattern. It then extracts a central portion
        of the result, effectively applying a padding factor to increase the resolution of the Fourier transform.

        Parameters
        ----------
        - scalar : numpy.ndarray
            A two-dimensional numpy array representing the scalar field of which the Fourier component
            is to be computed. This field could represent either the parallel or perpendicular projection of the far-field
            pattern onto the detector plane.

        Returns
        -------
        numpy.ndarray
            A two-dimensional numpy array representing the computed Fourier component. This array is a square
            section, extracted from the center of the full Fourier transform, with dimensions determined by the original
            sampling rate and the padding factor of the instance. The values in the array represent the intensity distribution
            of the light in the detector plane, providing insights into the spatial characteristics of the scattering pattern.

        The method uses numpy's fft module to perform the Fourier transform, applying a padding factor to the input to
        achieve a higher resolution in the Fourier domain. The resulting Fourier transform is then squared and fftshifted
        to center the zero-frequency component, and a central portion is extracted to match the intended output size.
        """
        # Calculate the target size based on the sampling and padding factor, and the indices for the central portion extraction.
        total_size = self.sampling * self.padding_factor
        offset = (total_size - self.sampling) // 2

        # Apply zero-padding to the scalar field to increase the resolution of the Fourier transform.
        padded_scalar = numpy.pad(scalar, pad_width=((offset, offset), (offset, offset)), mode='constant', constant_values=0)

        # Perform the two-dimensional inverse Fourier transform on the padded scalar field.
        fourier_transformed = numpy.fft.ifft2(padded_scalar)

        # Compute the squared magnitude and center the zero-frequency component.
        fourier_magnitude_squared = numpy.abs(numpy.fft.fftshift(fourier_transformed))**2

        # Extract the central portion corresponding to the original sampling rate adjusted by the padding factor.
        central_portion = fourier_magnitude_squared[offset:-offset, offset:-offset]

        return central_portion

    def plot(self, colormap: str = 'gray') -> None:
        """
        Plots the scatterer footprint using a 2D colormap.

        The method generates a plot representing the footprint of the scatterer, with the X and Y axes showing
        offset distances in micrometers, and the colormap depicting the mapping values.

        Parameters
        ----------
        colormap : str
            The colormap to use for the plot. Default is 'gray'.

        """
        with plt.style.context(MPSPlots.styles.mps):
            figure, ax = plt.subplots()

            ax.set(
                title='Scatterer Footprint',
                xlabel=r'Offset distance in X-axis [$\mu$m]',
                ylabel=r'Offset distance in Y-axis [$\mu$m]',
            )

            ax.pcolormesh(
                self.direct_y,
                self.direct_x,
                self.mapping,
                cmap=colormap
            )

            plt.show()
