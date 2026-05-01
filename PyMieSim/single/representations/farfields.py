#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Far-field visualization helpers for single-scatterer results."""

from typing import Sequence

import matplotlib.pyplot as pyplot
import numpy
from matplotlib import cm
from matplotlib import colors

from MPSPlots.colormaps import blue_black_red

from PyMieSim.units import ureg, Length


class FarFields:
    r"""
    Far-field scattering representation.

    This class computes and visualizes the complex far-field electric-field
    components returned by a single-scatterer setup.

    The far-field components are evaluated on a spherical angular mesh and are
    returned as the transverse components \(E_\phi\) and \(E_\theta\).

    The default visualization is a single three-dimensional spherical plot. A
    single field is shown at a time to keep rendering responsive and to make the
    figure easier to interpret.
    """

    def __init__(
        self,
        setup: object,
        sampling: int = 200,
        distance: Length = 1.0 * ureg.meter,
    ) -> None:
        self.setup = setup
        self.sampling = sampling
        self.distance = distance

        self.E_phi, self.E_theta, self.mesh = self.setup.get_farfields(
            sampling=self.sampling,
            distance=self.distance,
        )

    def plot(
        self,
        field: str = "total_intensity",
        unit_size: Sequence[float] = (5.0, 5.0),
        background_color: str = "white",
        colormap=None,
        opacity: float = 1.0,
        show_edges: bool = False,
        show_axis_label: bool = False,
        show_colorbar: bool = True,
        show_vector_field: bool = False,
        vector_count: int = 10,
        vector_radius: float = 1.1,
        vector_length: float = 0.08,
        vector_color: str = "black",
        scale: str = "log",
        percentile_clip: float | None = None,
        surface_linewidth: float = 0.0,
        surface_antialiased: bool = False,
        elevation: float = 25.0,
        azimuth: float = 35.0,
    ):
        r"""
        Plot one far-field quantity on a three-dimensional spherical surface.

        Parameters
        ----------
        field : str, optional
            Quantity to plot. Supported values are:

            ``"total_intensity"``
                Total intensity, computed as ``|E_phi|^2 + |E_theta|^2``.

            ``"phi_intensity"``
                Phi component intensity, computed as ``|E_phi|^2``.

            ``"theta_intensity"``
                Theta component intensity, computed as ``|E_theta|^2``.

            ``"phi_real"``
                Real part of \(E_\phi\).

            ``"phi_imag"``
                Imaginary part of \(E_\phi\).

            ``"theta_real"``
                Real part of \(E_\theta\).

            ``"theta_imag"``
                Imaginary part of \(E_\theta\).

            ``"phi_phase"``
                Phase of \(E_\phi\).

            ``"theta_phase"``
                Phase of \(E_\theta\).

            ``"phase_difference"``
                Wrapped phase difference ``arg(E_theta) - arg(E_phi)``.
        unit_size : Sequence[float], optional
            Figure size in inches, given as ``(width, height)``.
        background_color : str, optional
            Figure and axis background color.
        colormap : str or matplotlib.colors.Colormap, optional
            Colormap used for the selected field. If ``None``,
            ``MPSPlots.colormaps.blue_black_red`` is used.
        opacity : float, optional
            Surface opacity.
        show_edges : bool, optional
            If ``True``, draw surface mesh edges.
        show_axis_label : bool, optional
            If ``True``, display Cartesian axis labels, ticks, and panes.
        show_colorbar : bool, optional
            If ``True``, add a colorbar.
        show_vector_field : bool, optional
            If ``True``, overlay sparse local basis vectors. Phi fields receive
            phi vectors, theta fields receive theta vectors, and total fields do
            not receive a vector overlay.
        vector_count : int, optional
            Number of angular samples per direction used for vector overlays.
        vector_radius : float, optional
            Radius at which direction vectors are drawn.
        vector_length : float, optional
            Display length of direction vectors.
        vector_color : str, optional
            Color of the direction vectors.
        scale : str, optional
            Scaling for intensity fields. Supported values are ``"linear"`` and
            ``"log"``. Ignored for signed and phase fields.
        percentile_clip : float, optional
            If provided, upper color limits are clipped to this percentile.
        surface_linewidth : float, optional
            Line width used for surface edges.
        surface_antialiased : bool, optional
            Matplotlib antialiasing flag for surfaces.
        elevation : float, optional
            Matplotlib 3D view elevation angle in degrees.
        azimuth : float, optional
            Matplotlib 3D view azimuth angle in degrees.

        Returns
        -------
        matplotlib.figure.Figure
            Matplotlib figure containing the spherical far-field plot.
        """
        field = field.lower().strip()

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

        field_values, field_label, field_type, field_basis = self._get_field_data(field)

        figure = pyplot.figure(
            figsize=(unit_size[0], unit_size[1]),
            facecolor=background_color,
        )

        ax = figure.add_subplot(1, 1, 1, projection="3d")

        self._plot_field_surface(
            ax=ax,
            x_coordinates=x_coordinates,
            y_coordinates=y_coordinates,
            z_coordinates=z_coordinates,
            field=field_values,
            label=field_label,
            field_type=field_type,
            colormap=colormap,
            opacity=opacity,
            show_edges=show_edges,
            show_colorbar=show_colorbar,
            scale=scale,
            percentile_clip=percentile_clip,
            surface_linewidth=surface_linewidth,
            surface_antialiased=surface_antialiased,
        )

        if show_vector_field:
            if field_basis == "phi":
                self.add_phi_vector_to_3d_plot(
                    ax=ax,
                    n_points=vector_count,
                    radius=vector_radius,
                    length=vector_length,
                    color=vector_color,
                    opacity=1.0,
                )

            if field_basis == "theta":
                self.add_theta_vector_to_3d_plot(
                    ax=ax,
                    n_points=vector_count,
                    radius=vector_radius,
                    length=vector_length,
                    color=vector_color,
                    opacity=1.0,
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

    def plot_heatmap(
        self,
        kind: str = "intensity",
        unit_size: Sequence[float] = (4.2, 3.6),
        colormap=None,
        background_color: str = "white",
        show_colorbar: bool = True,
        scale: str = "log",
        percentile_clip: float | None = None,
        interpolation: str = "nearest",
    ):
        r"""
        Plot far-field quantities as two-dimensional angular heatmaps.

        Parameters
        ----------
        kind : str, optional
            Plot type. Supported values are ``"intensity"``, ``"complex"``,
            and ``"phase"``.
        unit_size : Sequence[float], optional
            Size of each panel in inches, given as ``(width, height)``.
        colormap : str or matplotlib.colors.Colormap, optional
            Colormap used for the plotted values. If ``None``,
            ``MPSPlots.colormaps.blue_black_red`` is used.
        background_color : str, optional
            Figure background color.
        show_colorbar : bool, optional
            If ``True``, add one colorbar to each panel.
        scale : str, optional
            Scaling for intensity plots. Supported values are ``"linear"`` and
            ``"log"``. Ignored for ``"complex"`` and ``"phase"``.
        percentile_clip : float, optional
            If provided, upper color limits are clipped to this percentile of
            the finite plotted values.
        interpolation : str, optional
            Interpolation mode passed to ``imshow``.

        Returns
        -------
        matplotlib.figure.Figure
            Matplotlib figure containing the far-field heatmaps.
        """
        kind = kind.lower().strip()

        if kind == "intensity":
            return self._plot_intensity_heatmaps(
                unit_size=unit_size,
                colormap=colormap,
                background_color=background_color,
                show_colorbar=show_colorbar,
                scale=scale,
                percentile_clip=percentile_clip,
                interpolation=interpolation,
            )

        if kind == "complex":
            return self._plot_complex_heatmaps(
                unit_size=unit_size,
                colormap=colormap,
                background_color=background_color,
                show_colorbar=show_colorbar,
                percentile_clip=percentile_clip,
                interpolation=interpolation,
            )

        if kind == "phase":
            return self._plot_phase_heatmaps(
                unit_size=unit_size,
                colormap=colormap,
                background_color=background_color,
                show_colorbar=show_colorbar,
                interpolation=interpolation,
            )

        raise ValueError(
            "Unknown far-field heatmap kind. Expected one of: "
            "'intensity', 'complex', or 'phase'."
        )

    def add_theta_vector_to_3d_plot(
        self,
        ax,
        n_points: int = 10,
        opacity: float = 1.0,
        radius: float = 1.0,
        length: float = 0.08,
        color: str = "black",
    ) -> None:
        r"""
        Add sparse unit vectors in the local \( \theta \) direction.

        This reproduces the previous PyVista call

        ``transform_vectors_sph_to_cart(theta=theta, phi=phi, r=radius, u=1, v=0, w=0)``

        where ``theta`` is the azimuthal angle in degrees and ``phi`` is the
        polar angle in degrees.
        """
        theta_mesh, phi_mesh = self._make_pyvista_style_angle_mesh(n_points)

        x_coordinates, y_coordinates, z_coordinates = self._spherical_to_cartesian(
            theta=theta_mesh,
            phi=phi_mesh,
            radius=radius,
        )

        theta_x = numpy.cos(phi_mesh) * numpy.cos(theta_mesh)
        theta_y = numpy.cos(phi_mesh) * numpy.sin(theta_mesh)
        theta_z = -numpy.sin(phi_mesh)

        ax.quiver(
            x_coordinates,
            y_coordinates,
            z_coordinates,
            theta_x,
            theta_y,
            theta_z,
            length=length,
            normalize=True,
            color=color,
            alpha=opacity,
            linewidth=0.8,
            arrow_length_ratio=0.4,
        )

    def add_phi_vector_to_3d_plot(
        self,
        ax,
        n_points: int = 10,
        opacity: float = 1.0,
        radius: float = 1.0,
        length: float = 0.08,
        color: str = "black",
    ) -> None:
        r"""
        Add sparse unit vectors in the local \( \phi \) direction.

        This reproduces the previous PyVista call

        ``transform_vectors_sph_to_cart(theta=theta, phi=phi, r=radius, u=0, v=1, w=0)``

        where ``theta`` is the azimuthal angle in degrees and ``phi`` is the
        polar angle in degrees.
        """
        theta_mesh, phi_mesh = self._make_pyvista_style_angle_mesh(n_points)

        x_coordinates, y_coordinates, z_coordinates = self._spherical_to_cartesian(
            theta=theta_mesh,
            phi=phi_mesh,
            radius=radius,
        )

        phi_x = -numpy.sin(theta_mesh)
        phi_y = numpy.cos(theta_mesh)
        phi_z = numpy.zeros_like(theta_mesh)

        ax.quiver(
            x_coordinates,
            y_coordinates,
            z_coordinates,
            phi_x,
            phi_y,
            phi_z,
            length=length,
            normalize=True,
            color=color,
            alpha=opacity,
            linewidth=0.8,
            arrow_length_ratio=0.4,
        )

    def _plot_field_surface(
        self,
        ax,
        x_coordinates: numpy.ndarray,
        y_coordinates: numpy.ndarray,
        z_coordinates: numpy.ndarray,
        field: numpy.ndarray,
        label: str,
        field_type: str,
        colormap,
        opacity: float,
        show_edges: bool,
        show_colorbar: bool,
        scale: str,
        percentile_clip: float | None,
        surface_linewidth: float,
        surface_antialiased: bool,
    ) -> None:
        """
        Plot one scalar field on the far-field spherical mesh.
        """
        field = numpy.asarray(field, dtype=float)

        colormap_object = self._resolve_colormap(
            colormap=colormap,
            field_type=field_type,
        )

        normalization = self._get_normalization(
            field=field,
            field_type=field_type,
            scale=scale,
            percentile_clip=percentile_clip,
        )

        facecolors = colormap_object(normalization(field))
        facecolors[..., -1] = opacity

        edgecolor = "black" if show_edges else "none"
        linewidth = surface_linewidth if show_edges else 0.0

        ax.plot_surface(
            x_coordinates,
            y_coordinates,
            z_coordinates,
            facecolors=facecolors,
            rstride=1,
            cstride=1,
            linewidth=linewidth,
            edgecolor=edgecolor,
            antialiased=surface_antialiased,
            shade=False,
        )

        ax.set_title(label)

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
                label=label,
            )

    def _plot_intensity_heatmaps(
        self,
        unit_size: Sequence[float],
        colormap,
        background_color: str,
        show_colorbar: bool,
        scale: str,
        percentile_clip: float | None,
        interpolation: str,
    ):
        """
        Plot intensity-like quantities on angular heatmaps.
        """
        e_phi = self._complex_field_array(self.E_phi)
        e_theta = self._complex_field_array(self.E_theta)

        phi_intensity = numpy.abs(e_phi) ** 2
        theta_intensity = numpy.abs(e_theta) ** 2
        total_intensity = phi_intensity + theta_intensity

        fields = [
            phi_intensity,
            theta_intensity,
            total_intensity,
        ]

        labels = [
            r"$|E_\phi|^2$",
            r"$|E_\theta|^2$",
            r"$|E_\phi|^2 + |E_\theta|^2$",
        ]

        return self._plot_heatmap_grid(
            fields=fields,
            labels=labels,
            field_type="intensity",
            unit_size=unit_size,
            colormap=colormap if colormap is not None else blue_black_red,
            background_color=background_color,
            show_colorbar=show_colorbar,
            scale=scale,
            percentile_clip=percentile_clip,
            interpolation=interpolation,
        )

    def _plot_complex_heatmaps(
        self,
        unit_size: Sequence[float],
        colormap,
        background_color: str,
        show_colorbar: bool,
        percentile_clip: float | None,
        interpolation: str,
    ):
        """
        Plot real and imaginary far-field components on angular heatmaps.
        """
        e_phi = self._complex_field_array(self.E_phi)
        e_theta = self._complex_field_array(self.E_theta)

        fields = [
            e_phi.real,
            e_phi.imag,
            e_theta.real,
            e_theta.imag,
        ]

        labels = [
            r"$\operatorname{Re}(E_\phi)$",
            r"$\operatorname{Im}(E_\phi)$",
            r"$\operatorname{Re}(E_\theta)$",
            r"$\operatorname{Im}(E_\theta)$",
        ]

        return self._plot_heatmap_grid(
            fields=fields,
            labels=labels,
            field_type="signed",
            unit_size=unit_size,
            colormap=colormap if colormap is not None else blue_black_red,
            background_color=background_color,
            show_colorbar=show_colorbar,
            scale="linear",
            percentile_clip=percentile_clip,
            interpolation=interpolation,
        )

    def _plot_phase_heatmaps(
        self,
        unit_size: Sequence[float],
        colormap,
        background_color: str,
        show_colorbar: bool,
        interpolation: str,
    ):
        """
        Plot far-field phases on angular heatmaps.
        """
        e_phi = self._complex_field_array(self.E_phi)
        e_theta = self._complex_field_array(self.E_theta)

        phi_phase = numpy.angle(e_phi)
        theta_phase = numpy.angle(e_theta)
        phase_difference = self._wrap_phase(theta_phase - phi_phase)

        fields = [
            phi_phase,
            theta_phase,
            phase_difference,
        ]

        labels = [
            r"$\arg(E_\phi)$",
            r"$\arg(E_\theta)$",
            r"$\arg(E_\theta) - \arg(E_\phi)$",
        ]

        return self._plot_heatmap_grid(
            fields=fields,
            labels=labels,
            field_type="phase",
            unit_size=unit_size,
            colormap=colormap if colormap is not None else blue_black_red,
            background_color=background_color,
            show_colorbar=show_colorbar,
            scale="linear",
            percentile_clip=None,
            interpolation=interpolation,
        )

    def _plot_heatmap_grid(
        self,
        fields: list[numpy.ndarray],
        labels: list[str],
        field_type: str,
        unit_size: Sequence[float],
        colormap,
        background_color: str,
        show_colorbar: bool,
        scale: str,
        percentile_clip: float | None,
        interpolation: str,
    ):
        """
        Plot a list of fields as angular heatmaps.
        """
        figure = pyplot.figure(
            figsize=(unit_size[0] * len(fields), unit_size[1]),
            facecolor=background_color,
        )

        axes = [
            figure.add_subplot(1, len(fields), index + 1)
            for index in range(len(fields))
        ]

        for ax, field, label in zip(axes, fields, labels):
            image = self._plot_single_heatmap(
                ax=ax,
                field=field,
                label=label,
                field_type=field_type,
                colormap=colormap,
                scale=scale,
                percentile_clip=percentile_clip,
                interpolation=interpolation,
                extent=(0.0, 360.0, 180.0, 0.0),
            )

            ax.set_xlabel(r"$\theta$ [deg]")
            ax.set_ylabel(r"$\phi$ [deg]")
            ax.set_facecolor(background_color)

            if show_colorbar:
                figure.colorbar(
                    image,
                    ax=ax,
                    fraction=0.046,
                    pad=0.04,
                    label=label,
                )

        figure.tight_layout()
        pyplot.show()

        return figure

    def _plot_single_heatmap(
        self,
        ax,
        field: numpy.ndarray,
        label: str,
        field_type: str,
        colormap,
        scale: str,
        percentile_clip: float | None,
        interpolation: str,
        extent: tuple[float, float, float, float],
    ):
        """
        Plot a single angular heatmap.
        """
        field = numpy.asarray(field, dtype=float)

        colormap_object = self._resolve_colormap(
            colormap=colormap,
            field_type=field_type,
        )

        normalization = self._get_normalization(
            field=field,
            field_type=field_type,
            scale=scale,
            percentile_clip=percentile_clip,
        )

        image = ax.imshow(
            field,
            origin="upper",
            extent=extent,
            aspect="auto",
            interpolation=interpolation,
            cmap=colormap_object,
            norm=normalization,
        )

        ax.set_title(label)

        return image

    def _get_field_data(
        self,
        field: str,
    ) -> tuple[numpy.ndarray, str, str, str | None]:
        """
        Return data, label, field type, and basis family for one named field.
        """
        e_phi = self._complex_field_array(self.E_phi)
        e_theta = self._complex_field_array(self.E_theta)

        phi_intensity = numpy.abs(e_phi) ** 2
        theta_intensity = numpy.abs(e_theta) ** 2
        total_intensity = phi_intensity + theta_intensity

        field_map = {
            "total_intensity": (
                total_intensity,
                r"$|E_\phi|^2 + |E_\theta|^2$",
                "intensity",
                None,
            ),
            "phi_intensity": (
                phi_intensity,
                r"$|E_\phi|^2$",
                "intensity",
                "phi",
            ),
            "theta_intensity": (
                theta_intensity,
                r"$|E_\theta|^2$",
                "intensity",
                "theta",
            ),
            "phi_real": (
                e_phi.real,
                r"$\operatorname{Re}(E_\phi)$",
                "signed",
                "phi",
            ),
            "phi_imag": (
                e_phi.imag,
                r"$\operatorname{Im}(E_\phi)$",
                "signed",
                "phi",
            ),
            "theta_real": (
                e_theta.real,
                r"$\operatorname{Re}(E_\theta)$",
                "signed",
                "theta",
            ),
            "theta_imag": (
                e_theta.imag,
                r"$\operatorname{Im}(E_\theta)$",
                "signed",
                "theta",
            ),
            "phi_phase": (
                numpy.angle(e_phi),
                r"$\arg(E_\phi)$",
                "phase",
                "phi",
            ),
            "theta_phase": (
                numpy.angle(e_theta),
                r"$\arg(E_\theta)$",
                "phase",
                "theta",
            ),
            "phase_difference": (
                self._wrap_phase(numpy.angle(e_theta) - numpy.angle(e_phi)),
                r"$\arg(E_\theta) - \arg(E_\phi)$",
                "phase",
                None,
            ),
        }

        if field not in field_map:
            raise ValueError(
                "Unknown far-field field. Expected one of: "
                f"{', '.join(field_map)}."
            )

        return field_map[field]

    def _get_normalization(
        self,
        field: numpy.ndarray,
        field_type: str,
        scale: str,
        percentile_clip: float | None,
    ):
        """
        Build the appropriate Matplotlib normalization for a field.
        """
        finite_values = field[numpy.isfinite(field)]

        if finite_values.size == 0:
            return colors.Normalize(vmin=0.0, vmax=1.0)

        if field_type == "phase":
            return colors.Normalize(vmin=-numpy.pi, vmax=numpy.pi)

        if field_type == "signed":
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

        if field_type == "intensity":
            positive_values = finite_values[finite_values > 0.0]

            if positive_values.size == 0:
                return colors.Normalize(vmin=0.0, vmax=1.0)

            upper_limit = numpy.nanmax(positive_values)

            if percentile_clip is not None:
                upper_limit = numpy.nanpercentile(
                    positive_values,
                    percentile_clip,
                )

            if not numpy.isfinite(upper_limit) or upper_limit <= 0.0:
                upper_limit = 1.0

            if scale == "log":
                lower_limit = numpy.nanmax(positive_values) * 1e-6
                lower_limit = max(lower_limit, numpy.nanmin(positive_values))

                if lower_limit >= upper_limit:
                    lower_limit = upper_limit * 1e-6

                return colors.LogNorm(
                    vmin=lower_limit,
                    vmax=upper_limit,
                )

            if scale == "linear":
                return colors.Normalize(
                    vmin=0.0,
                    vmax=upper_limit,
                )

            raise ValueError(
                "Invalid intensity scale. Expected 'linear' or 'log'."
            )

        return colors.Normalize(
            vmin=numpy.nanmin(finite_values),
            vmax=numpy.nanmax(finite_values),
        )

    def _resolve_colormap(
        self,
        colormap,
        field_type: str,
    ):
        """
        Return a Matplotlib colormap object.
        """
        if colormap is None:
            colormap = blue_black_red

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
        Apply common formatting to one 3D axis.
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

    def _set_equal_axis_limits(self, ax) -> None:
        """
        Set symmetric equal limits on a Matplotlib 3D axis.
        """
        axis_limit = 1.15

        ax.set_xlim(-axis_limit, axis_limit)
        ax.set_ylim(-axis_limit, axis_limit)
        ax.set_zlim(-axis_limit, axis_limit)

        if hasattr(ax, "set_box_aspect"):
            ax.set_box_aspect((1.0, 1.0, 1.0))

    def _make_pyvista_style_angle_mesh(
        self,
        n_points: int,
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """
        Build the same angular sampling convention used by the previous PyVista
        vector overlay.

        PyVista calls used:

        ``theta = numpy.linspace(0, 360, n_points)``
        ``phi = numpy.linspace(180, 0, n_points)``

        The arrays returned here are in radians.
        """
        theta_degrees = numpy.linspace(0.0, 360.0, n_points)
        phi_degrees = numpy.linspace(180.0, 0.0, n_points)

        theta_radians = numpy.deg2rad(theta_degrees)
        phi_radians = numpy.deg2rad(phi_degrees)

        return numpy.meshgrid(
            theta_radians,
            phi_radians,
            indexing="ij",
        )

    def _spherical_to_cartesian(
        self,
        theta: numpy.ndarray,
        phi: numpy.ndarray,
        radius: float,
    ) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
        """
        Convert PyVista-style spherical coordinates to Cartesian coordinates.

        ``theta`` is azimuthal and ``phi`` is polar.
        """
        x_coordinates = radius * numpy.sin(phi) * numpy.cos(theta)
        y_coordinates = radius * numpy.sin(phi) * numpy.sin(theta)
        z_coordinates = radius * numpy.cos(phi)

        return x_coordinates, y_coordinates, z_coordinates

    def _complex_field_array(
        self,
        value,
    ) -> numpy.ndarray:
        """
        Convert a complex Pint quantity or array-like object to a square array.

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
        value: numpy.ndarray,
    ) -> numpy.ndarray:
        """
        Return a two-dimensional square array compatible with angular plotting.

        The backend may return structured quantities either as ``(N, N)`` arrays
        or as flattened arrays with ``N * N`` values. Flattened values are
        reshaped using Fortran ordering to match the previous PyVista flattening
        convention.
        """
        array = numpy.asarray(value)

        if array.ndim == 2:
            return array

        flat_array = array.ravel()
        expected_size = self.sampling * self.sampling

        if flat_array.size != expected_size:
            raise ValueError(
                "Cannot reshape array to the structured far-field mesh. "
                f"Expected {expected_size} values from sampling={self.sampling}, "
                f"but received {flat_array.size}."
            )

        return flat_array.reshape(
            (self.sampling, self.sampling),
            order="F",
        )

    def _wrap_phase(
        self,
        phase: numpy.ndarray,
    ) -> numpy.ndarray:
        """
        Wrap a phase-difference array to the interval ``[-pi, pi]``.

        ``numpy.angle`` already returns individual phases in this interval. This
        helper is used only for phase differences, where direct subtraction can
        otherwise create artificial jumps near the branch cut.
        """
        return (phase + numpy.pi) % (2.0 * numpy.pi) - numpy.pi
