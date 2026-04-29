#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Sequence

import matplotlib.pyplot as pyplot
import numpy
from matplotlib import cm
from matplotlib import colors

from PyMieSim.units import ureg
from PyMieSim.utils import spherical_to_cartesian
from PyMieSim.mesh import FullMesh  # Necessary for loading the class, even if not directly used in this file


class SPF:
    r"""
    Scattering phase function representation.

    The scattering phase function describes how the scattered intensity varies
    as a function of the scattering direction around the particle.

    The SPF is sampled on the backend spherical mesh and can be rendered either
    as a normalized radial surface or as an intensity map on a unit sphere.

    Parameters
    ----------
    setup : object
        Single scatterer setup object exposing ``get_spf``.
    sampling : int, optional
        Number of angular samples used in each spherical mesh direction.

    Attributes
    ----------
    setup : object
        Simulation setup used to compute the scattering phase function.
    sampling : int
        Angular sampling used for the structured SPF mesh.
    SPF : numpy.ndarray
        Scattering phase function sampled on the spherical mesh.
    mesh : object
        Full spherical mesh returned by the backend.
    """

    def __init__(
        self,
        setup,
        sampling: int = 200,
    ) -> None:
        self.setup = setup
        self.sampling = sampling

        self.SPF, self.mesh = self.setup.get_spf(
            sampling=self.sampling,
            distance=1.0 * ureg.meter,
        )

    def plot(
        self,
        unit_size: Sequence[float] = (5.0, 5.0),
        background_color: str = "white",
        show_edges: bool = False,
        edge_color: str = "black",
        colormap="viridis",
        opacity: float = 1.0,
        set_surface: bool = True,
        show_axis_label: bool = False,
        scale: str = "linear",
        percentile_clip: float | None = None,
        surface_linewidth: float = 0.1,
        surface_antialiased: bool = False,
        elevation: float = 25.0,
        azimuth: float = 35.0,
    ):
        r"""
        Plot the scattering phase function on a Matplotlib 3D surface.

        Parameters
        ----------
        unit_size : Sequence[float], optional
            Figure size in inches, given as ``(width, height)``.
        background_color : str, optional
            Figure and axis background color.
        show_edges : bool, optional
            If ``True``, draw surface mesh edges. Edges are also drawn
            automatically when ``surface_linewidth`` is greater than zero.
        edge_color : str, optional
            Color used for surface mesh edges.
        colormap : str or matplotlib.colors.Colormap, optional
            Colormap used for the SPF intensity values.
        opacity : float, optional
            Surface opacity.
        set_surface : bool, optional
            If ``True``, the surface radius is proportional to the normalized
            SPF intensity. If ``False``, the SPF intensity is drawn on a unit sphere.
        show_axis_label : bool, optional
            If ``True``, display Cartesian axis labels, ticks, and panes.
            If ``False``, hide the axis frame.
        scale : str, optional
            Color scaling. Supported values are ``"linear"`` and ``"log"``.
        percentile_clip : float, optional
            If provided, the upper color limit is clipped to this percentile of
            the finite positive values.
        surface_linewidth : float, optional
            Line width used for surface edges. Values greater than zero enable
            mesh edges even if ``show_edges`` is ``False``.
        surface_antialiased : bool, optional
            Matplotlib antialiasing flag for surfaces.
        elevation : float, optional
            Matplotlib 3D view elevation angle in degrees.
        azimuth : float, optional
            Matplotlib 3D view azimuth angle in degrees.

        Returns
        -------
        matplotlib.figure.Figure
            Matplotlib figure containing the SPF visualization.
        """
        figure = pyplot.figure(
            figsize=(unit_size[0], unit_size[1]),
            facecolor=background_color,
        )

        ax = figure.add_subplot(1, 1, 1, projection="3d")

        self._add_to_3d_ax(
            ax=ax,
            set_surface=set_surface,
            show_edges=show_edges,
            edge_color=edge_color,
            colormap=colormap,
            opacity=opacity,
            scale=scale,
            percentile_clip=percentile_clip,
            surface_linewidth=surface_linewidth,
            surface_antialiased=surface_antialiased,
        )

        self._format_3d_axis(
            ax=ax,
            background_color=background_color,
            show_axis_label=show_axis_label,
            elevation=elevation,
            azimuth=azimuth,
        )

        figure.tight_layout()
        pyplot.show()

        return figure

    def _add_to_3d_ax(
        self,
        ax,
        set_surface: bool = True,
        show_edges: bool = False,
        edge_color: str = "black",
        colormap="viridis",
        opacity: float = 1.0,
        scale: str = "linear",
        percentile_clip: float | None = None,
        surface_linewidth: float = 0.0,
        surface_antialiased: bool = False,
    ) -> None:
        r"""
        Add the SPF surface to a Matplotlib 3D axis.

        The surface geometry follows the scattering phase function if
        ``set_surface`` is ``True``. In all cases, the surface colors are linked
        to the SPF intensity itself.
        """
        intensity = self._spf_array()

        if set_surface:
            maximum_intensity = numpy.nanmax(intensity)

            if not numpy.isfinite(maximum_intensity) or maximum_intensity <= 0.0:
                radius = numpy.ones_like(intensity)
            else:
                radius = intensity / maximum_intensity
        else:
            radius = numpy.ones_like(intensity)

        x_coordinates, y_coordinates, z_coordinates = spherical_to_cartesian(
            r=radius,
            phi=self.mesh.spherical_mesh.phi.to("radian").magnitude,
            theta=self.mesh.spherical_mesh.theta.to("radian").magnitude,
        )

        x_coordinates = self._as_square_array(x_coordinates)
        y_coordinates = self._as_square_array(y_coordinates)
        z_coordinates = self._as_square_array(z_coordinates)

        intensity = self._as_square_array(intensity)

        colormap_object = self._resolve_colormap(colormap)

        normalization = self._get_normalization(
            intensity=intensity,
            scale=scale,
            percentile_clip=percentile_clip,
        )

        facecolors = colormap_object(normalization(intensity))
        facecolors[..., -1] = opacity

        draw_edges = show_edges or surface_linewidth > 0.0
        effective_edge_color = edge_color if draw_edges else "none"
        effective_linewidth = surface_linewidth if draw_edges else 0.0

        ax.plot_surface(
            x_coordinates,
            y_coordinates,
            z_coordinates,
            facecolors=facecolors,
            rstride=1,
            cstride=1,
            linewidth=effective_linewidth,
            edgecolor=effective_edge_color,
            antialiased=surface_antialiased,
            shade=False,
        )

    def _get_normalization(
        self,
        intensity: numpy.ndarray,
        scale: str,
        percentile_clip: float | None,
    ):
        """
        Build the Matplotlib normalization for the SPF intensity values.
        """
        finite_values = intensity[numpy.isfinite(intensity)]

        if finite_values.size == 0:
            return colors.Normalize(vmin=0.0, vmax=1.0)

        non_negative_values = finite_values[finite_values >= 0.0]

        if non_negative_values.size == 0:
            return colors.Normalize(vmin=0.0, vmax=1.0)

        upper_limit = numpy.nanmax(non_negative_values)

        if percentile_clip is not None:
            upper_limit = numpy.nanpercentile(
                non_negative_values,
                percentile_clip,
            )

        if not numpy.isfinite(upper_limit) or upper_limit <= 0.0:
            upper_limit = 1.0

        if scale == "linear":
            return colors.Normalize(
                vmin=0.0,
                vmax=upper_limit,
            )

        if scale == "log":
            strictly_positive_values = non_negative_values[non_negative_values > 0.0]

            if strictly_positive_values.size == 0:
                return colors.Normalize(vmin=0.0, vmax=upper_limit)

            lower_limit = numpy.nanmax(strictly_positive_values) * 1e-6
            lower_limit = max(lower_limit, numpy.nanmin(strictly_positive_values))

            if lower_limit >= upper_limit:
                lower_limit = upper_limit * 1e-6

            return colors.LogNorm(
                vmin=lower_limit,
                vmax=upper_limit,
            )

        raise ValueError(
            "Invalid SPF scale. Expected 'linear' or 'log'."
        )

    def _resolve_colormap(
        self,
        colormap,
    ):
        """
        Return a Matplotlib colormap object.
        """
        if isinstance(colormap, str):
            return cm.get_cmap(colormap)

        return colormap

    def _format_3d_axis(
        self,
        ax,
        background_color: str,
        show_axis_label: bool,
        elevation: float,
        azimuth: float,
    ) -> None:
        """
        Apply common formatting to the 3D axis.
        """
        ax.set_facecolor(background_color)
        ax.view_init(elev=elevation, azim=azimuth)

        self._set_equal_axis_limits(ax)

        if show_axis_label:
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
        else:
            ax.set_axis_off()

    def _set_equal_axis_limits(
        self,
        ax,
    ) -> None:
        """
        Set symmetric equal limits on a Matplotlib 3D axis.
        """
        axis_limit = 1.15

        ax.set_xlim(-axis_limit, axis_limit)
        ax.set_ylim(-axis_limit, axis_limit)
        ax.set_zlim(-axis_limit, axis_limit)

        if hasattr(ax, "set_box_aspect"):
            ax.set_box_aspect((1.0, 1.0, 1.0))

    def _spf_array(
        self,
    ) -> numpy.ndarray:
        """
        Return the SPF as a square NumPy array.
        """
        return self._as_square_array(
            self._quantity_to_magnitude_array(self.SPF)
        )

    def _quantity_to_magnitude_array(
        self,
        value,
    ) -> numpy.ndarray:
        """
        Convert a Pint quantity or array like object to a NumPy array.
        """
        if hasattr(value, "magnitude"):
            return numpy.asarray(value.magnitude)

        return numpy.asarray(value)

    def _as_square_array(
        self,
        value,
    ) -> numpy.ndarray:
        """
        Return a two dimensional square array compatible with surface plotting.
        """
        array = numpy.asarray(value)

        if array.ndim == 2:
            return array

        flat_array = array.ravel()
        expected_size = self.sampling * self.sampling

        if flat_array.size != expected_size:
            raise ValueError(
                "Cannot reshape array to the structured SPF mesh. "
                f"Expected {expected_size} values from sampling={self.sampling}, "
                f"but received {flat_array.size}."
            )

        return flat_array.reshape(
            (self.sampling, self.sampling),
            order="F",
        )