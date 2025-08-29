#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Tuple
import numpy
import MPSPlots
from pydantic.dataclasses import dataclass
import matplotlib.pyplot as plt
from TypedUnit import Length

from PyMieSim.utils import config_dict


@dataclass(config=config_dict, kw_only=True)
class NearField:
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
    x_range : Tuple[Length, Length]
        Range of x coordinates (x_min, x_max) in meters.
    y_range : Tuple[Length, Length]
        Range of y coordinates (y_min, y_max) in meters.
    z : Tuple[Length, Length]
        Range of z coordinates (z_min, z_max) in meters, or single z value for 2D slice.
    resolution : int or tuple
        Number of points along each axis. If int, same resolution for all axes.
    field_components : list
        List of field components to compute: ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "\|E\|", "\|H\|"]

    Raises
    ------
    RuntimeError
        If scatterer doesn't support near-field computation (e.g., cylinders).
    ValueError
        If coordinate ranges or field components are invalid.
    """

    scatterer: object
    x_range: Tuple[Length, Length]
    y_range: Tuple[Length, Length]
    z: Length
    sampling: int
    field_components: list[str]

    def __post_init__(self):
        # Set up coordinate grids
        self.radius = self.scatterer.diameter / 2
        self._setup_coordinates()

        # Compute all requested field components
        self._compute_fields()

    def _setup_coordinates(self):
        """Set up coordinate grids for field computation."""
        # Handle resolution
        self.nx = self.sampling
        self.ny = self.sampling

        # Create coordinate arrays
        self.x = numpy.linspace(self.x_range[0], self.x_range[1], self.nx)
        self.y = numpy.linspace(self.y_range[0], self.y_range[1], self.ny)

        # Create coordinate meshes
        self.X, self.Y = numpy.meshgrid(self.x, self.y, indexing="ij")
        self.Z = numpy.full_like(self.X, self.z)

    def _compute_fields(self):
        """Compute all requested field components."""
        # Compute field using C++ backend
        self.fields = {}
        for component in self.field_components:
            self.fields[component] = self.scatterer._cpp_compute_nearfields(
                x=self.X.flatten().to("meter").magnitude,
                y=self.Y.flatten().to("meter").magnitude,
                z=self.Z.flatten().to("meter").magnitude,
                field_type=component,
            ).reshape(self.X.shape)

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

    def plot(self, colormap: str = "viridis", figure_size: tuple = (6, 6)) -> None:
        """
        Plot 2D near-field distribution.

        Parameters
        ----------
        colormap : str
            Matplotlib colormap name.
        figure_size : tuple
            Figure size (width, height).
        """
        with plt.style.context(MPSPlots.styles.mps):
            fig, axes = plt.subplots(
                figsize=figure_size, nrows=1, ncols=len(self.fields), squeeze=False
            )

        for component, ax in zip(self.fields.keys(), axes.flatten()):
            field_data = numpy.abs(self.fields[component]).astype(float)

            print(self.x, self.x.__class__)
            # Plot field
            im = ax.pcolormesh(
                self.x.magnitude,
                self.y.magnitude,
                field_data,
                cmap=colormap,
            )

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label(f"{component} field")

            # Show scatterer boundary
            circle = plt.Circle(
                (0, 0), self.radius.magnitude, fill=False, color="white", linewidth=2
            )
            ax.add_patch(circle)

            ax.set_title(f"Near-field {component} distribution")
            ax.set_aspect("equal")

        plt.grid(False)
        # plt.tight_layout()
        plt.show()

    def plot_enhancement(
        self,
        field_type: str = "electric",
        enhancement_threshold: float = 2.0,
        colormap: str = "plasma",
        figsize: tuple = (10, 8),
    ) -> None:
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

            enhancement_data = enhancement[:, :, 0]
            extent = [
                self.x[0] * 1e9,
                self.x[-1] * 1e9,
                self.y[0] * 1e9,
                self.y[-1] * 1e9,
            ]

            # Plot enhancement
            im = ax.imshow(
                enhancement_data.T,
                extent=extent,
                origin="lower",
                cmap=colormap,
                aspect="equal",
            )

            # Highlight hot-spots
            hot_spots = enhancement_data > enhancement_threshold
            if numpy.any(hot_spots):
                ax.contour(
                    enhancement_data.T,
                    levels=[enhancement_threshold],
                    colors="red",
                    linewidths=2,
                    extent=extent,
                )

            # Show scatterer boundary
            circle = plt.Circle(
                (0, 0), self.radius, fill=False, color="black", linewidth=2
            )
            ax.add_patch(circle)

            # Colorbar and labels
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label(f"{field_type.title()} field enhancement")

            ax.set_xlabel("x [nm]")
            ax.set_ylabel("y [nm]")
            ax.set_title(f"{field_type.title()} field enhancement")

            plt.tight_layout()
            plt.show()
