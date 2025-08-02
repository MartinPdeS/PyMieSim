#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import MPSPlots
from pydantic.dataclasses import dataclass
from PyMieSim.special_functions import rotate_on_x
import matplotlib.pyplot as plt
from PyMieSim.units import meter
from PyMieSim.single.representations.base import BaseRepresentation, config_dict


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


@dataclass(config=config_dict, kw_only=True)
class NearField():
    r"""
    Represents electromagnetic near-field distributions around a scatterer.

    This class provides comprehensive near-field visualization capabilities for spherical
    and core-shell scatterers using the internal (cn, dn) and external (an, bn)
    multipole coefficients from Mie theory.

    The electromagnetic fields are computed using vector spherical harmonics:

    .. math::
        \mathbf{E}(r) = \begin{cases}
            \sum_{n=1}^{N} [c_n \mathbf{M}_{n}^{(1)} + d_n \mathbf{N}_{n}^{(1)}] & r < a \\
            \sum_{n=1}^{N} [a_n \mathbf{M}_{n}^{(3)} + b_n \mathbf{N}_{n}^{(3)}] & r > a
        \end{cases}

    Where :math:`\mathbf{M}_n` and :math:`\mathbf{N}_n` are vector spherical harmonics,
    and superscripts (1) and (3) denote regular and outgoing wave solutions.

    Parameters
    ----------
    scatterer : BaseScatterer
        The scatterer object with computed cn/dn coefficients.
    x_range : tuple
        Range of x coordinates (x_min, x_max) in meters.
    y_range : tuple
        Range of y coordinates (y_min, y_max) in meters.
    z_range : tuple or float
        Range of z coordinates (z_min, z_max) in meters, or single z value for 2D slice.
    resolution : int or tuple
        Number of points along each axis. If int, same resolution for all axes.
    field_components : list
        List of field components to compute: ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "|E|", "|H|"]

    Raises
    ------
    RuntimeError
        If scatterer doesn't support near-field computation (e.g., cylinders).
    ValueError
        If coordinate ranges or field components are invalid.

    Examples
    --------
    >>> import numpy as np
    >>> from PyMieSim.single.scatterer import Sphere
    >>> from PyMieSim.single.source import PlaneWave
    >>> from PyMieSim.single.representations import NearField
    >>> from PyMieSim.units import nanometer, RIU
    >>>
    >>> # Create sphere with near-field capability
    >>> source = PlaneWave(wavelength=532*nanometer, polarization=0)
    >>> sphere = Sphere(diameter=200*nanometer, property=1.5*RIU, source=source)
    >>> sphere.compute_cn_dn()  # Enable near-field computation
    >>>
    >>> # Compute near-field in a 2D slice
    >>> near_field = NearField(
    ...     scatterer=sphere,
    ...     x_range=(-200e-9, 200e-9),
    ...     y_range=(-200e-9, 200e-9),
    ...     z_range=0,  # 2D slice at z=0
    ...     resolution=100,
    ...     field_components=["Ex", "Ey", "|E|"]
    ... )
    >>>
    >>> # Visualize the results
    >>> near_field.plot_2d()
    >>> near_field.plot_enhancement()
    """
    scatterer: object
    x_range: tuple
    y_range: tuple
    z_range: tuple | float
    resolution: int | tuple = 100
    field_components: list = None

    def __post_init__(self):
        # Set default field components
        if self.field_components is None:
            self.field_components = ["Ex", "Ey", "Ez", "|E|"]

        # Validate scatterer compatibility
        if not hasattr(self.scatterer, 'compute_near_field_py'):
            raise RuntimeError(
                "Scatterer does not support near-field computation. "
                "Currently supported: Sphere, CoreShell with computed cn/dn coefficients."
            )

        # Determine if this is 2D or 3D computation
        self.is_2d = isinstance(self.z_range, (int, float))

        # Set up coordinate grids
        self._setup_coordinates()

        # Compute all requested field components
        self.fields = {}
        self._compute_fields()

        # Determine scatterer radius for boundary visualization
        self.radius = getattr(self.scatterer, 'radius',
                             getattr(self.scatterer, 'diameter', 0) / 2)

    def _setup_coordinates(self):
        """Set up coordinate grids for field computation."""
        # Handle resolution
        if isinstance(self.resolution, int):
            self.nx = self.ny = self.nz = self.resolution
        else:
            self.nx, self.ny, self.nz = self.resolution

        # Create coordinate arrays
        self.x = numpy.linspace(self.x_range[0], self.x_range[1], self.nx)
        self.y = numpy.linspace(self.y_range[0], self.y_range[1], self.ny)

        if self.is_2d:
            self.z = numpy.array([self.z_range])
            self.nz = 1
        else:
            self.z = numpy.linspace(self.z_range[0], self.z_range[1], self.nz)

        # Create coordinate meshes
        if self.is_2d:
            self.X, self.Y = numpy.meshgrid(self.x, self.y, indexing='ij')
            self.Z = numpy.full_like(self.X, self.z_range)
        else:
            self.X, self.Y, self.Z = numpy.meshgrid(self.x, self.y, self.z, indexing='ij')

    def _compute_fields(self):
        """Compute all requested field components."""
        for component in self.field_components:
            if component not in ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "|E|", "|H|"]:
                raise ValueError(f"Invalid field component: {component}")

            # Flatten coordinates for C++ computation
            x_flat = self.X.flatten()
            y_flat = self.Y.flatten()
            z_flat = self.Z.flatten()

            # Compute field using C++ backend
            field_values = self.scatterer.compute_near_field_py(
                x_flat, y_flat, z_flat, component, self.radius
            )

            # Reshape to grid
            self.fields[component] = field_values.reshape(self.X.shape)

    def get_field_enhancement(self, field_type="electric"):
        """
        Compute field enhancement factor |E_total|/|E_incident|.

        Parameters
        ----------
        field_type : str
            "electric" or "magnetic" field enhancement.

        Returns
        -------
        numpy.ndarray
            Field enhancement map.
        """
        if field_type == "electric":
            if "|E|" not in self.fields:
                raise ValueError("Electric field magnitude |E| not computed")
            return numpy.abs(self.fields["|E|"])
        elif field_type == "magnetic":
            if "|H|" not in self.fields:
                raise ValueError("Magnetic field magnitude |H| not computed")
            return numpy.abs(self.fields["|H|"])
        else:
            raise ValueError("field_type must be 'electric' or 'magnetic'")

    def plot_2d(self,
                component: str = "|E|",
                plane: str = "xy",
                slice_position: float = 0,
                colormap: str = 'viridis',
                show_scatterer: bool = True,
                enhancement_levels: list = None,
                figsize: tuple = (10, 8)) -> None:
        """
        Plot 2D near-field distribution.

        Parameters
        ----------
        component : str
            Field component to plot.
        plane : str
            Plane to plot: "xy", "xz", or "yz".
        slice_position : float
            Position of slice through 3D data (ignored for 2D data).
        colormap : str
            Matplotlib colormap name.
        show_scatterer : bool
            Whether to show scatterer boundary.
        enhancement_levels : list
            Contour levels for field enhancement.
        figsize : tuple
            Figure size (width, height).
        """
        if component not in self.fields:
            raise ValueError(f"Component {component} not computed")

        with plt.style.context(MPSPlots.styles.mps):
            fig, ax = plt.subplots(figsize=figsize)

            # Get field data for plotting
            if self.is_2d and plane == "xy":
                field_data = numpy.abs(self.fields[component][:, :, 0])
                extent = [self.x[0]*1e9, self.x[-1]*1e9, self.y[0]*1e9, self.y[-1]*1e9]
                xlabel, ylabel = "x [nm]", "y [nm]"
            else:
                # Handle 3D slicing
                field_data = self._get_slice_data(component, plane, slice_position)
                extent = self._get_extent(plane)
                xlabel, ylabel = self._get_labels(plane)

            # Plot field
            im = ax.imshow(
                field_data.T,
                extent=extent,
                origin='lower',
                cmap=colormap,
                aspect='equal'
            )

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label(f'{component} field')

            # Show scatterer boundary
            if show_scatterer and plane == "xy":
                circle = plt.Circle((0, 0), self.radius*1e9,
                                  fill=False, color='white', linewidth=2)
                ax.add_patch(circle)

            # Add enhancement contours
            if enhancement_levels and "|E|" in self.fields:
                enhancement = self.get_field_enhancement()
                if self.is_2d:
                    enhancement_slice = enhancement[:, :, 0]
                else:
                    enhancement_slice = self._get_slice_data("|E|", plane, slice_position)

                ax.contour(extent[:2], extent[2:], enhancement_slice.T,
                          levels=enhancement_levels, colors='white', alpha=0.7)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(f'Near-field {component} distribution')

            plt.tight_layout()
            plt.show()

    def plot_enhancement(self,
                        field_type: str = "electric",
                        enhancement_threshold: float = 2.0,
                        colormap: str = 'plasma',
                        figsize: tuple = (10, 8)) -> None:
        """
        Plot field enhancement map with hot-spots highlighted.

        Parameters
        ----------
        field_type : str
            "electric" or "magnetic" field enhancement.
        enhancement_threshold : float
            Threshold for highlighting hot-spots.
        colormap : str
            Matplotlib colormap name.
        figsize : tuple
            Figure size (width, height).
        """
        enhancement = self.get_field_enhancement(field_type)

        with plt.style.context(MPSPlots.styles.mps):
            fig, ax = plt.subplots(figsize=figsize)

            if self.is_2d:
                enhancement_data = enhancement[:, :, 0]
                extent = [self.x[0]*1e9, self.x[-1]*1e9, self.y[0]*1e9, self.y[-1]*1e9]
            else:
                # Take central slice for 3D data
                mid_z = self.nz // 2
                enhancement_data = enhancement[:, :, mid_z]
                extent = [self.x[0]*1e9, self.x[-1]*1e9, self.y[0]*1e9, self.y[-1]*1e9]

            # Plot enhancement
            im = ax.imshow(
                enhancement_data.T,
                extent=extent,
                origin='lower',
                cmap=colormap,
                aspect='equal'
            )

            # Highlight hot-spots
            hot_spots = enhancement_data > enhancement_threshold
            if numpy.any(hot_spots):
                ax.contour(enhancement_data.T, levels=[enhancement_threshold],
                          colors='red', linewidths=2,
                          extent=extent)

            # Show scatterer boundary
            circle = plt.Circle((0, 0), self.radius*1e9,
                              fill=False, color='black', linewidth=2)
            ax.add_patch(circle)

            # Colorbar and labels
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label(f'{field_type.title()} field enhancement')

            ax.set_xlabel('x [nm]')
            ax.set_ylabel('y [nm]')
            ax.set_title(f'{field_type.title()} field enhancement')

            plt.tight_layout()
            plt.show()

    def plot_3d(self,
                component: str = "|E|",
                threshold: float = None,
                colormap: str = 'viridis',
                opacity: float = 0.8,
                show_scatterer: bool = True) -> None:
        """
        Plot 3D near-field isosurfaces using PyVista.

        Parameters
        ----------
        component : str
            Field component to plot.
        threshold : float
            Isosurface threshold value. If None, uses 50% of maximum.
        colormap : str
            PyVista colormap name.
        opacity : float
            Isosurface opacity.
        show_scatterer : bool
            Whether to show scatterer boundary.
        """
        if self.is_2d:
            raise ValueError("3D plotting requires 3D field data")

        if component not in self.fields:
            raise ValueError(f"Component {component} not computed")

        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("PyVista required for 3D plotting: pip install pyvista")

        # Create structured grid
        grid = pv.StructuredGrid(self.X, self.Y, self.Z)
        field_data = numpy.abs(self.fields[component])
        grid.point_data[component] = field_data.flatten(order='F')

        # Set threshold
        if threshold is None:
            threshold = 0.5 * field_data.max()

        # Create plotter
        plotter = pv.Plotter()

        # Add isosurface
        isosurface = grid.contour([threshold])
        plotter.add_mesh(isosurface, cmap=colormap, opacity=opacity)

        # Add scatterer boundary
        if show_scatterer:
            sphere_boundary = pv.Sphere(radius=self.radius, center=(0, 0, 0))
            plotter.add_mesh(sphere_boundary, color='gray', opacity=0.3)

        # Customize plot
        plotter.add_axes()
        plotter.add_scalar_bar(title=f'{component} field')
        plotter.set_background('white')

        plotter.show()

    def _get_slice_data(self, component: str, plane: str, position: float):
        """Extract 2D slice from 3D field data."""
        field_data = numpy.abs(self.fields[component])

        if plane == "xy":
            # Find closest z index
            z_idx = numpy.argmin(numpy.abs(self.z - position))
            return field_data[:, :, z_idx]
        elif plane == "xz":
            y_idx = numpy.argmin(numpy.abs(self.y - position))
            return field_data[:, y_idx, :]
        elif plane == "yz":
            x_idx = numpy.argmin(numpy.abs(self.x - position))
            return field_data[x_idx, :, :]
        else:
            raise ValueError("plane must be 'xy', 'xz', or 'yz'")

    def _get_extent(self, plane: str):
        """Get plot extent for given plane."""
        if plane == "xy":
            return [self.x[0]*1e9, self.x[-1]*1e9, self.y[0]*1e9, self.y[-1]*1e9]
        elif plane == "xz":
            return [self.x[0]*1e9, self.x[-1]*1e9, self.z[0]*1e9, self.z[-1]*1e9]
        elif plane == "yz":
            return [self.y[0]*1e9, self.y[-1]*1e9, self.z[0]*1e9, self.z[-1]*1e9]

    def _get_labels(self, plane: str):
        """Get axis labels for given plane."""
        if plane == "xy":
            return "x [nm]", "y [nm]"
        elif plane == "xz":
            return "x [nm]", "z [nm]"
        elif plane == "yz":
            return "y [nm]", "z [nm]"
