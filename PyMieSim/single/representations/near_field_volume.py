#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import MPSPlots
from pydantic.dataclasses import dataclass
import matplotlib.pyplot as plt
from PyMieSim.single.representations.base import config_dict
from PyMieSim import units

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
    x_range: tuple[units.Quantity, units.Quantity]
    y_range: tuple[units.Quantity, units.Quantity]
    z_range: tuple[units.Quantity, units.Quantity] | units.Quantity
    resolution: int | tuple = 100
    field_components: list = None

    def __post_init__(self):
        # Set default field components
        if self.field_components is None:
            self.field_components = ["Ex", "Ey", "Ez", "|E|"]

        # Determine if this is 2D or 3D computation
        self.is_2d = isinstance(self.z_range, (int, float))

        # Set up coordinate grids
        self.radius = self.scatterer.diameter / 2
        self._setup_coordinates()

        # Compute all requested field components
        self.fields = {}
        self._compute_fields()

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
            field_values = self.scatterer._cpp_compute_near_field(
                x_flat, y_flat, z_flat, component, self.radius.to('meter').magnitude
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
                extent = [self.x[0] * 1e9, self.x[-1] * 1e9, self.y[0] * 1e9, self.y[-1]*1e9]
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

            # # Show scatterer boundary
            # if show_scatterer and plane == "xy":
            #     circle = plt.Circle((0, 0), self.radius * 1e9, fill=False, color='white', linewidth=2)
            #     ax.add_patch(circle)

            # Add enhancement contours
            if enhancement_levels and "|E|" in self.fields:
                enhancement = self.get_field_enhancement()
                if self.is_2d:
                    enhancement_slice = enhancement[:, :, 0]
                else:
                    enhancement_slice = self._get_slice_data("|E|", plane, slice_position)

                ax.contour(extent[:2], extent[2:], enhancement_slice.T, levels=enhancement_levels, colors='white', alpha=0.7)

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
