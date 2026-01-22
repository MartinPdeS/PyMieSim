#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from typing import Tuple, Optional, Sequence, Dict, List
import numpy
import matplotlib.pyplot as plt
from MPSPlots.colormaps import blue_black_red
import MPSPlots

from PyMieSim.units import Length


def _normalize_vector(vector: numpy.ndarray, eps: float = 1e-15) -> numpy.ndarray:
    """
    Normalize a 3D vector.

    Parameters
    ----------
    vector
        Vector to normalize.
    eps
        Threshold below which the norm is considered zero.

    Returns
    -------
    numpy.ndarray
        Normalized vector.
    """
    norm = float(numpy.linalg.norm(vector))
    if norm < eps:
        raise ValueError("Vector norm too small to normalize")
    return vector / norm


def _make_plane_basis(plane_normal: Sequence[float]) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """
    Build an orthonormal basis for a plane.

    Parameters
    ----------
    plane_normal
        Normal vector of the plane in 3D.

    Returns
    -------
    Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        (n_hat, u_hat, v_hat) where u_hat and v_hat span the plane.
    """
    n_hat = _normalize_vector(numpy.array(plane_normal, dtype=float))

    if abs(n_hat[2]) < 0.9:
        candidate = numpy.array([0.0, 0.0, 1.0], dtype=float)
    else:
        candidate = numpy.array([1.0, 0.0, 0.0], dtype=float)

    u_hat = numpy.cross(n_hat, candidate)
    u_hat = _normalize_vector(u_hat)

    v_hat = numpy.cross(n_hat, u_hat)
    v_hat = _normalize_vector(v_hat)

    return n_hat, u_hat, v_hat


class NearField:
    """
    Near field visualization for a scatterer.

    This class samples fields over an arbitrary plane defined by an origin and a normal.
    All coordinates are kept as unit aware objects because the C++ wrapper accepts units.

    Important concept
    -----------------
    The plot is parameterized by two in plane coordinates (u, v):

        r(u, v) = r0 + u * u_hat + v * v_hat

    The ranges u_range and v_range define the sampled extent in that plane.
    These replace the old idea of x_range and y_range, because the plane is not necessarily a z slice.

    Parameters
    ----------
    scatterer
        Object exposing compute_total_nearfields, compute_incident_nearfields, compute_scattered_nearfields.
    """

    def __init__(self, scatterer: object) -> None:
        """
        Initialize NearField object.

        Parameters
        ----------
        scatterer
            Object exposing compute_total_nearfields, compute_incident_nearfields, compute_scattered_nearfields.
        """
        self.scatterer = scatterer

        self.sampling: int = 200

        self.u = None
        self.v = None
        self.U = None
        self.V = None

        self.X = None
        self.Y = None
        self.Z = None

        self.plane_origin = None
        self.plane_normal = None
        self.plane_u_hat = None
        self.plane_v_hat = None

        self.u_range = None
        self.v_range = None

    def _default_uv_range(self, extent_scale: float) -> Tuple[Tuple[Length, Length], Tuple[Length, Length]]:
        """
        Define a default sampling extent in the plane.

        The default is symmetric about zero and based on the particle radius.
        extent_scale lets you zoom out.

        Parameters
        ----------
        extent_scale
            Multiplies the default half width.

        Returns
        -------
        Tuple[Tuple[Length, Length], Tuple[Length, Length]]
            (u_range, v_range)
        """
        if not hasattr(self.scatterer, "diameter"):
            raise ValueError("scatterer must expose diameter to infer a default plotting extent")

        radius = 0.5 * self.scatterer.diameter
        half_width = extent_scale * radius

        u_range = (-half_width, half_width)
        v_range = (-half_width, half_width)

        return u_range, v_range

    def _setup_coordinates_plane(
        self,
        *,
        sampling: int,
        plane_origin: Optional[Sequence[Length]] = None,
        plane_normal: Optional[Sequence[float]] = None,
        u_range: Optional[Tuple[Length, Length]] = None,
        v_range: Optional[Tuple[Length, Length]] = None,
        extent_scale: float = 2.5,
    ) -> None:
        """
        Setup the coordinate grid for the sampling plane.

        Parameters
        ----------
        sampling
            Number of points along each in plane direction (u and v).
        plane_origin
            Origin of the sampling plane in 3D space. If None, defaults to (0, 0, z).
        plane_normal
            Normal vector of the sampling plane. If None, defaults to (0, 0, 1).
        u_range
            Range of plane coordinate u along u_hat. If None, inferred from scatterer diameter and extent_scale.
        v_range
            Range of plane coordinate v along v_hat. If None, inferred from scatterer diameter and extent_scale.
        extent_scale
            Scaling factor applied to the default extent derived from the scatterer radius.
        """
        self.sampling = sampling

        if u_range is None or v_range is None:
            inferred_u_range, inferred_v_range = self._default_uv_range(extent_scale=extent_scale)
            if u_range is None:
                u_range = inferred_u_range
            if v_range is None:
                v_range = inferred_v_range

        self.u_range = u_range
        self.v_range = v_range

        u_min, u_max = u_range
        v_min, v_max = v_range

        self.u = numpy.linspace(u_min, u_max, sampling)
        self.v = numpy.linspace(v_min, v_max, sampling)

        U, V = numpy.meshgrid(self.u, self.v, indexing="ij")
        self.U = U
        self.V = V

        origin = (plane_origin[0], plane_origin[1], plane_origin[2])

        if plane_normal is None:
            normal = (0.0, 0.0, 1.0)
        else:
            normal = (float(plane_normal[0]), float(plane_normal[1]), float(plane_normal[2]))

        _, u_hat, v_hat = _make_plane_basis(normal)

        self.plane_origin = origin
        self.plane_normal = normal
        self.plane_u_hat = u_hat
        self.plane_v_hat = v_hat

        origin_x, origin_y, origin_z = origin

        self.X = origin_x + U * u_hat[0] + V * v_hat[0]
        self.Y = origin_y + U * u_hat[1] + V * v_hat[1]
        self.Z = origin_z + U * u_hat[2] + V * v_hat[2]

    def _compute_fields(self, field_components: List[str], *, type: str) -> Dict[str, numpy.ndarray]:
        """
        Compute near fields over the sampling plane.

        Parameters
        ----------
        field_components
            List of field components to compute, for example Ex, Ey, Ez, |E|.
        type
            Type of field to compute: total, incident, or scattered.

        Returns
        -------
        Dict[str, numpy.ndarray]
            Dictionary mapping component name to 2D array of values on the plane grid.
        """
        fields: Dict[str, numpy.ndarray] = {}

        x_flat = self.X.reshape(-1)
        y_flat = self.Y.reshape(-1)
        z_flat = self.Z.reshape(-1)

        for component in field_components:
            print(f"Computing {type} near field component: {component}")

            if type == "total":
                values = self.scatterer.compute_total_nearfields(x=x_flat, y=y_flat, z=z_flat, field_type=component)
            elif type == "incident":
                values = self.scatterer.compute_incident_nearfields(x=x_flat, y=y_flat, z=z_flat, field_type=component)
            elif type == "scattered":
                values = self.scatterer.compute_scattered_nearfields(x=x_flat, y=y_flat, z=z_flat, field_type=component)
            else:
                raise ValueError(f"Unknown field type: {type}")

            fields[component] = numpy.asarray(values).reshape(self.X.shape)

        return fields

    @staticmethod
    def _parse_field_components(field_components: Sequence[str]) -> Tuple[List[str], List[str]]:
        """
        Parse field components into their main component and part.

        Parameters
        ----------
        field_components
            Sequence of field components, for example Ex:real, |E|:abs.

        Returns
        -------
        Tuple[List[str], List[str]]
            Two lists: one with the main components and one with the parts.
        """
        components: List[str] = []
        parts: List[str] = []
        for item in field_components:
            if ":" in item:
                component, part = item.split(":")
            else:
                component, part = item, "real"
            components.append(component)
            parts.append(part)
        return components, parts

    def _add_scatterer_outline_on_plane(self, ax) -> None:
        """
        Draw the intersection of the sphere with the current plane.

        If the plane intersects the sphere, the intersection is a circle in plane coordinates (u, v).
        The circle center in (u, v) is the projection of the sphere center onto the plane.
        """
        if self.U is None or self.V is None:
            return

        if not hasattr(self.scatterer, "diameter"):
            return

        radius = 0.5 * self.scatterer.diameter

        sphere_center = numpy.array([0.0, 0.0, 0.0], dtype=float)

        origin_x, origin_y, origin_z = self.plane_origin
        plane_origin = numpy.array([float(origin_x), float(origin_y), float(origin_z)], dtype=float)

        n_hat = _normalize_vector(numpy.array(self.plane_normal, dtype=float))
        u_hat = numpy.array(self.plane_u_hat, dtype=float)
        v_hat = numpy.array(self.plane_v_hat, dtype=float)

        signed_distance = float(numpy.dot((sphere_center - plane_origin), n_hat))
        abs_distance = abs(signed_distance)

        if abs_distance > radius:
            return

        intersection_radius = numpy.sqrt(max(radius * radius - abs_distance * abs_distance, 0.0))

        closest_point = sphere_center - signed_distance * n_hat
        relative = closest_point - plane_origin

        u0 = float(numpy.dot(relative, u_hat))
        v0 = float(numpy.dot(relative, v_hat))

        circle = plt.Circle(
            (u0, v0),
            intersection_radius,
            fill=False,
            color="white",
            linewidth=2.0,
        )
        ax.add_patch(circle)

    @MPSPlots.helper.post_mpl_plot
    def plot(
        self,
        *field_components: str,
        type: str = "scattered",
        sampling: int = 200,
        plane_origin: Optional[Sequence[Length]] = None,
        plane_normal: Optional[Sequence[float]] = None,
        u_range: Optional[Tuple[Length, Length]] = None,
        v_range: Optional[Tuple[Length, Length]] = None,
        extent_scale: float = 2.5,
        colormap: str = blue_black_red,
        figure_size: tuple = (6, 6),
        show_scatterer_outline: bool = True,
    ):
        """
        Plot near fields over a specified plane.

        Parameters
        ----------
        field_components
            List of field components to plot, for example Ex:real, |E|:abs.
        type
            Type of field to compute: total, incident, or scattered.
        sampling
            Number of points along each in plane direction (u and v).
        plane_origin
            Origin of the sampling plane in 3D space. If None, defaults to (0, 0, z).
        plane_normal
            Normal vector of the sampling plane. If None, defaults to (0, 0, 1).
        u_range
            Range of plane coordinate u along u_hat. If None, inferred from scatterer diameter and extent_scale.
        v_range
            Range of plane coordinate v along v_hat. If None, inferred from scatterer diameter and extent_scale.
        extent_scale
            Scaling factor applied to the default extent derived from the scatterer radius.
        colormap
            Colormap to use for plotting.
        figure_size
            Size of the matplotlib figure.
        show_scatterer_outline
            If True, draw the intersection circle of the sphere with the plane.
        """
        if len(field_components) == 0:
            raise ValueError("Provide at least one field component, for example Ex:real or |E|:abs")

        self._setup_coordinates_plane(
            sampling=sampling,
            plane_origin=plane_origin,
            plane_normal=plane_normal,
            u_range=u_range,
            v_range=v_range,
            extent_scale=extent_scale,
        )

        components, parts = self._parse_field_components(field_components)
        fields = self._compute_fields(list(components), type=type)

        figure, axes = plt.subplots(
            figsize=figure_size,
            nrows=1,
            ncols=len(fields),
            squeeze=False,
        )

        for component, part, ax in zip(fields.keys(), parts, axes.flatten()):
            values = fields[component]

            if part == "real":
                field_data = values.real.astype(float)
            elif part == "imag":
                field_data = values.imag.astype(float)
            elif part == "abs":
                field_data = numpy.abs(values).astype(float)
            else:
                raise ValueError(f"Unknown field component part: {part}")

            im = ax.pcolormesh(
                self.u,
                self.v,
                field_data,
                cmap=colormap,
                shading="auto",
            )

            if show_scatterer_outline:
                self._add_scatterer_outline_on_plane(ax)

            plt.colorbar(im, ax=ax)

            ax.set_aspect("equal")
            ax.set_xlabel("plane u")
            ax.set_ylabel("plane v")
            ax.set_title(f"{type} {component} {part}")
            ax.grid(visible=False, which="both", axis="both")

        return figure
