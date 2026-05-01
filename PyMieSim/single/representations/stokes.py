#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Stokes-parameter visualization helpers for single-scatterer results."""

from typing import Sequence

import matplotlib.pyplot as pyplot
import numpy
from matplotlib import cm
from matplotlib import colors

from MPSPlots.colormaps import blue_black_red

from PyMieSim.units import ureg
from PyMieSim.mesh import FullMesh  # Necessary for loading the class, even if not directly used in this file


class Stokes:
    r"""
    Stokes parameter representation.

    The Stokes parameters describe the polarization state of the scattered light.

    .. math::

        I &= |E_x|^2 + |E_y|^2 \\
        Q &= |E_x|^2 - |E_y|^2 \\
        U &= 2 \operatorname{Re}\{E_x E_y^*\} \\
        V &= 2 \operatorname{Im}\{E_x E_y^*\}

    Parameters
    ----------
    setup : object
        Single-scatterer setup object exposing ``get_stokes``.
    sampling : int, optional
        Number of angular samples used in each spherical mesh direction.
    distance : pint.Quantity, optional
        Observation distance at which the Stokes parameters are evaluated.

    Attributes
    ----------
    setup : object
        Simulation setup used to compute the Stokes parameters.
    sampling : int
        Angular sampling used for the structured Stokes mesh.
    distance : pint.Quantity
        Observation distance.
    I : pint.Quantity
        Total scattered intensity.
    Q : pint.Quantity
        Linear polarization contrast between the x and y axes.
    U : pint.Quantity
        Linear polarization contrast between the +45 and -45 degree axes.
    V : pint.Quantity
        Circular polarization contrast.
    mesh : object
        Full spherical mesh returned by the backend.
    """

    def __init__(
        self,
        setup: object,
        sampling: int = 200,
        distance: ureg.Quantity = 1.0 * ureg.meter,
    ) -> None:
        self.setup = setup
        self.sampling = sampling
        self.distance = distance

        self.I, self.Q, self.U, self.V, self.mesh = self.setup.get_stokes(
            sampling=self.sampling,
            distance=self.distance,
        )

    def plot(
        self,
        field: str = "I",
        unit_size: Sequence[float] = (5.0, 5.0),
        background_color: str = "white",
        show_edges: bool = False,
        edge_color: str = "black",
        colormap=blue_black_red,
        opacity: float = 1.0,
        show_axis_label: bool = False,
        show_colorbar: bool = True,
        percentile_clip: float | None = None,
        surface_linewidth: float = 0.0,
        surface_antialiased: bool = False,
        elevation: float = 25.0,
        azimuth: float = 35.0,
    ):
        r"""
        Plot one Stokes parameter on a Matplotlib 3D spherical surface.

        Examples
        --------
        .. code-block:: python

            stokes.plot("I")
            stokes.plot("Q")
            stokes.plot("U")
            stokes.plot("V")

        Parameters
        ----------
        field : str, optional
            Stokes parameter to plot. Supported values are ``"I"``, ``"Q"``,
            ``"U"``, and ``"V"``.
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
            Colormap used for the Stokes parameter values.
        opacity : float, optional
            Surface opacity.
        show_axis_label : bool, optional
            If ``True``, display Cartesian axis labels, ticks, and panes.
            If ``False``, hide the axis frame.
        show_colorbar : bool, optional
            If ``True``, add a colorbar.
        percentile_clip : float, optional
            If provided, symmetric color limits are clipped to this percentile
            of the finite absolute values.
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
            Matplotlib figure containing the selected Stokes parameter.
        """
        field_values, field_label = self._get_field_data(field)

        cartesian = self.mesh.spherical_mesh.to_cartesian()

        x_coordinates = self._as_square_array(
            self._quantity_to_magnitude_array(cartesian.x, "meter")
        )
        y_coordinates = self._as_square_array(
            self._quantity_to_magnitude_array(cartesian.y, "meter")
        )
        z_coordinates = self._as_square_array(
            self._quantity_to_magnitude_array(cartesian.z, "meter")
        )

        figure = pyplot.figure(
            figsize=(unit_size[0], unit_size[1]),
            facecolor=background_color,
        )

        ax = figure.add_subplot(1, 1, 1, projection="3d")

        self._add_to_3d_ax(
            ax=ax,
            x_coordinates=x_coordinates,
            y_coordinates=y_coordinates,
            z_coordinates=z_coordinates,
            field=field_values,
            field_label=field_label,
            show_edges=show_edges,
            edge_color=edge_color,
            colormap=colormap,
            opacity=opacity,
            show_colorbar=show_colorbar,
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
        x_coordinates: numpy.ndarray,
        y_coordinates: numpy.ndarray,
        z_coordinates: numpy.ndarray,
        field: numpy.ndarray,
        field_label: str,
        show_edges: bool = False,
        edge_color: str = "black",
        colormap=blue_black_red,
        opacity: float = 1.0,
        show_colorbar: bool = True,
        percentile_clip: float | None = None,
        surface_linewidth: float = 0.0,
        surface_antialiased: bool = False,
    ) -> None:
        r"""
        Add one Stokes field to a Matplotlib 3D axis.

        The scalar ordering follows the previous PyVista implementation, where
        each field was attached to the ``StructuredGrid`` as
        ``field.flatten(order="F")``.
        """
        colormap_object = self._resolve_colormap(colormap)

        normalization = self._get_normalization(
            field=field,
            percentile_clip=percentile_clip,
        )

        facecolors = colormap_object(normalization(field))
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

        ax.set_title(f"{field_label} field")

        if show_colorbar:
            scalar_mappable = cm.ScalarMappable(
                norm=normalization,
                cmap=colormap_object,
            )
            scalar_mappable.set_array([])

            ax.figure.colorbar(
                scalar_mappable,
                ax=ax,
                fraction=0.046,
                pad=0.04,
                label=f"{field_label} field",
            )

    def _get_field_data(
        self,
        field: str,
    ) -> tuple[numpy.ndarray, str]:
        """
        Return the selected Stokes parameter and its display label.
        """
        field_key = field.upper().strip()

        field_map = {
            "I": self.I,
            "Q": self.Q,
            "U": self.U,
            "V": self.V,
        }

        if field_key not in field_map:
            raise ValueError(
                "Unknown Stokes field. Expected one of: 'I', 'Q', 'U', or 'V'."
            )

        return self._stokes_array(field_map[field_key]), field_key

    def _get_normalization(
        self,
        field: numpy.ndarray,
        percentile_clip: float | None,
    ):
        """
        Build a Matplotlib normalization for a Stokes parameter.

        ``I`` is non-negative in normal use, while ``Q``, ``U``, and ``V`` can be
        signed. A symmetric normalization around zero is therefore used for all
        fields so that the color semantics remain consistent across parameters.
        """
        finite_values = field[numpy.isfinite(field)]

        if finite_values.size == 0:
            return colors.Normalize(vmin=-1.0, vmax=1.0)

        maximum_absolute_value = numpy.nanmax(numpy.abs(finite_values))

        if percentile_clip is not None:
            maximum_absolute_value = numpy.nanpercentile(
                numpy.abs(finite_values),
                percentile_clip,
            )

        if not numpy.isfinite(maximum_absolute_value) or maximum_absolute_value <= 0.0:
            maximum_absolute_value = 1.0

        return colors.Normalize(
            vmin=-maximum_absolute_value,
            vmax=maximum_absolute_value,
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

    def _stokes_array(
        self,
        value,
    ) -> numpy.ndarray:
        """
        Return a Stokes parameter as a square NumPy array.

        The transpose reproduces the scalar ordering used by the previous
        PyVista implementation, where fields were attached to the
        ``StructuredGrid`` as ``field.flatten(order="F")``.
        """
        return self._as_square_array(
            self._quantity_to_magnitude_array(value)
        ).T

    def _quantity_to_magnitude_array(
        self,
        value,
        unit: str | None = None,
    ) -> numpy.ndarray:
        """
        Convert a Pint quantity or array-like object to a NumPy array.
        """
        if hasattr(value, "to") and unit is not None:
            return numpy.asarray(value.to(unit).magnitude)

        if hasattr(value, "magnitude"):
            return numpy.asarray(value.magnitude)

        return numpy.asarray(value)

    def _as_square_array(
        self,
        value,
    ) -> numpy.ndarray:
        """
        Return a two-dimensional square array compatible with surface plotting.
        """
        array = numpy.asarray(value)

        if array.ndim == 2:
            return array

        flat_array = array.ravel()
        expected_size = self.sampling * self.sampling

        if flat_array.size != expected_size:
            raise ValueError(
                "Cannot reshape array to the structured Stokes mesh. "
                f"Expected {expected_size} values from sampling={self.sampling}, "
                f"but received {flat_array.size}."
            )

        return flat_array.reshape(
            (self.sampling, self.sampling),
            order="F",
        )
