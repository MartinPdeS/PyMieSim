from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Optional, Sequence, Tuple, Type, Union

import numpy as np
import pyvista as pv

from MPSPlots.colormaps import blue_black_red

from PyMieSim.single.scatterer import Sphere, CoreShell, Cylinder
from PyMieSim.single.source import Gaussian, PlaneWave
from PyMieSim.single.detector import Photodiode, CoherentMode


PyVistaColor = Union[str, Tuple[float, float, float], Sequence[float]]


@dataclass
class SystemPlotter:
    """
    3D visualization helper for PyMieSim objects using PyVista.

    The class provides a minimal, extensible registry that maps object types
    (e.g., `Sphere`, `Gaussian`, `Photodiode`) to plotting callbacks that add
    corresponding meshes to a shared `pyvista.Plotter`.

    Parameters
    ----------
    colormap : Any, optional
        Colormap used for scalar fields (e.g., coherent collecting field).
        Defaults to `blue_black_red`.
    show_axis_label : bool, optional
        If True, show axis labels at the origin. Defaults to True.
    context_sphere_radius : float, optional
        Radius of a translucent context sphere. If <= 0, no context sphere
        is added. Defaults to 1.0.
    context_sphere_opacity : float, optional
        Opacity of the context sphere. Defaults to 0.3.

    Notes
    -----
    This plotter intentionally uses simple, schematic geometry for most objects
    (sphere, cylinder, arrow, cone) to provide a quick visual overview of the
    system layout rather than an exact physical rendering.
    """
    show_axis_label: bool = True
    context_sphere_radius: float = 1.0
    context_sphere_opacity: float = 0.3

    plotter_factory: Callable[[], pv.Plotter] = pv.Plotter

    _registry: Dict[Type[Any], Callable[[Any, pv.Plotter], None]] = field(
        default_factory=dict, init=False, repr=False
    )

    def __post_init__(self) -> None:
        self.colormap = blue_black_red
        # Order matters: more specific types first when you have inheritance.
        self.register(Sphere, self._add_sphere)
        self.register(CoreShell, self._add_sphere)
        self.register(Cylinder, self._add_cylinder)

        self.register(Gaussian, self._add_gaussian)
        self.register(PlaneWave, self._add_planewave)

        self.register(CoherentMode, self._add_coherent_detector)
        self.register(Photodiode, self._add_photodiode)

    def register(self, object_type: Type[Any], handler: Callable[[Any, pv.Plotter], None]) -> None:
        """
        Register a plotting handler for a given object type.

        Parameters
        ----------
        object_type : type
            Type (class) that the handler can plot.
        handler : callable
            Function with signature `(obj, scene) -> None` that adds meshes to `scene`.
        """
        self._registry[object_type] = handler

    def plot(
        self,
        *objects: Any,
        data: Optional[Any] = None,
        scene: Optional[pv.Plotter] = None,
        add_context_sphere: bool = True,
        show: bool = True,
    ) -> pv.Plotter:
        """
        Plot a system of objects in a shared 3D scene.

        Parameters
        ----------
        *objects : Any
            Objects to be plotted (e.g., scatterers, sources, detectors).
            Supported types are registered internally; extend via `register`.
        data : Any, optional
            Optional object that provides `_add_to_3d_ax(scene=..., colormap=...)`.
            If provided, it will be added after system objects.
        scene : pyvista.Plotter, optional
            Existing scene to add meshes to. If None, a new plotter is created.
        add_context_sphere : bool, optional
            If True, adds a translucent sphere of radius `context_sphere_radius`.
            Defaults to True.
        show : bool, optional
            If True, calls `scene.show()`. Defaults to True.

        Returns
        -------
        pyvista.Plotter
            The plotter used for the rendering (useful if `show=False`).

        Raises
        ------
        TypeError
            If an object type is not registered.
        AttributeError
            If `data` is provided but lacks `_add_to_3d_ax`.
        """
        active_scene = scene if scene is not None else self.plotter_factory()

        for obj in objects:
            handler = self._resolve_handler(obj)
            if handler is None:
                raise TypeError(
                    f"Unsupported object type: {type(obj).__name__}. "
                    f"Register a handler via `SystemPlotter.register(...)`."
                )
            handler(obj, active_scene)

        if data is not None:
            self._add_data(data=data, scene=active_scene)

        if add_context_sphere and self.context_sphere_radius > 0.0:
            sphere = pv.Sphere(radius=float(self.context_sphere_radius))
            active_scene.add_mesh(sphere, opacity=float(self.context_sphere_opacity))

        active_scene.add_axes_at_origin(
            zlabel="",
            xlabel="",
            ylabel="",
            labels_off=not self.show_axis_label,
        )

        if show:
            active_scene.show()

        return active_scene

    def _resolve_handler(self, obj: Any) -> Optional[Callable[[Any, pv.Plotter], None]]:
        # Resolve based on MRO so subclasses are handled without explicit registration.
        for cls in type(obj).mro():
            handler = self._registry.get(cls)
            if handler is not None:
                return handler
        return None

    @staticmethod
    def _add_sphere(obj: Any, scene: pv.Plotter, color: PyVistaColor = "black", opacity: float = 1.0) -> None:
        """
        Add a schematic sphere scatterer.

        Parameters
        ----------
        obj : Any
            Sphere like scatterer object (unused, kept for signature uniformity).
        scene : pyvista.Plotter
            Target plotter.
        color : str or tuple, optional
            Mesh color. Defaults to "black".
        opacity : float, optional
            Mesh opacity. Defaults to 1.0.
        """
        shape = pv.Sphere(center=(0.0, 0.0, 0.0), radius=0.1, theta_resolution=100, phi_resolution=100)
        scene.add_mesh(shape, color=color, opacity=opacity)

    @staticmethod
    def _add_cylinder(obj: Any, scene: pv.Plotter, color: PyVistaColor = "black", opacity: float = 1.0) -> None:
        """
        Add a schematic infinite cylinder scatterer.

        Parameters
        ----------
        obj : Any
            Cylinder scatterer object (unused, kept for signature uniformity).
        scene : pyvista.Plotter
            Target plotter.
        color : str or tuple, optional
            Mesh color. Defaults to "black".
        opacity : float, optional
            Mesh opacity. Defaults to 1.0.
        """
        shape = pv.Cylinder(
            center=(0.0, 0.0, 0.0),
            radius=0.1,
            height=2.0,
            direction=(0.0, -1.0, 0.0),
            resolution=100,
        )
        scene.add_mesh(shape, color=color, opacity=opacity)

    @staticmethod
    def _numerical_aperture_to_half_angle_radian(numerical_aperture: float) -> float:
        """
        Convert numerical aperture to a half angle (radians) for visualization.

        Parameters
        ----------
        numerical_aperture : float
            Numerical aperture magnitude.

        Returns
        -------
        float
            Half angle in radians.

        Notes
        -----
        This mirrors the behavior in your original function, including the
        fallback when NA > 1 (often indicating an effective NA in a medium).
        """
        if numerical_aperture <= 1.0:
            return float(np.arcsin(numerical_aperture))
        return float(np.arcsin(numerical_aperture - 1.0) + np.pi / 2.0)

    def _add_gaussian(self, obj: Gaussian, scene: pv.Plotter, color: PyVistaColor = "red") -> None:
        """
        Add a schematic Gaussian beam cone using the source NA.

        Parameters
        ----------
        obj : PyMieSim.single.source.Gaussian
            Gaussian source instance. Must provide `NA` with a `.magnitude`.
        scene : pyvista.Plotter
            Target plotter.
        color : str or tuple, optional
            Mesh color. Defaults to "red".
        """
        numerical_aperture = float(obj.NA.magnitude)
        half_angle_rad = self._numerical_aperture_to_half_angle_radian(numerical_aperture)
        half_angle_deg = float(np.rad2deg(half_angle_rad))

        cone_mesh = pv.Cone(
            center=(0.0, 0.0, -0.5),
            direction=(0.0, 0.0, 1.0),
            height=0.9,
            resolution=100,
            angle=half_angle_deg,
        )
        scene.add_mesh(cone_mesh, color=color, opacity=0.3)

    @staticmethod
    def _add_planewave(obj: PlaneWave, scene: pv.Plotter, color: PyVistaColor = "blue", opacity: float = 1.0) -> None:
        """
        Add a schematic plane wave direction arrow.

        Parameters
        ----------
        obj : PyMieSim.single.source.PlaneWave
            PlaneWave source instance (unused, kept for signature uniformity).
        scene : pyvista.Plotter
            Target plotter.
        color : str or tuple, optional
            Mesh color. Defaults to "blue".
        opacity : float, optional
            Mesh opacity. Defaults to 1.0.
        """
        shape = pv.Arrow(
            start=(0.0, 0.0, -1.0),
            direction=(0.0, 0.0, 1.0),
            tip_length=0.2,
            tip_radius=0.05,
            shaft_radius=0.02,
        )
        scene.add_mesh(shape, color=color, opacity=opacity)

    def _add_coherent_detector(
        self,
        obj: CoherentMode,
        scene: pv.Plotter,
        cone_color: PyVistaColor = "blue",
        field_point_size: float = 20.0,
    ) -> None:
        """
        Visualize a coherent mode detector collecting field.

        Parameters
        ----------
        obj : PyMieSim.single.detector.CoherentMode
            Coherent detector instance providing `_cpp_mesh` (with cartesian x/y/z)
            and `scalar_field` (complex scalar field values).
        scene : pyvista.Plotter
            Target plotter.
        cone_color : str or tuple, optional
            Color of the collection cone. Defaults to "blue".
        field_point_size : float, optional
            Point size for the field visualization. Defaults to 20.0.
        """
        coordinates = np.vstack(
            (
                obj._cpp_mesh.cartesian.x,
                obj._cpp_mesh.cartesian.y,
                obj._cpp_mesh.cartesian.z,
            )
        )

        points = pv.wrap(coordinates.T)

        scalar_field_real = np.asarray(obj.scalar_field).real
        abs_max = float(np.abs(scalar_field_real).max()) if scalar_field_real.size else 1.0

        mapping = scene.add_points(
            points,
            scalars=scalar_field_real,
            point_size=float(field_point_size),
            render_points_as_spheres=True,
            cmap=self.colormap,
            show_scalar_bar=False,
            clim=[-abs_max, abs_max],
        )

        mean_vector = coordinates.mean(axis=1)
        cone_mesh = pv.Cone(
            center=mean_vector / 2.0,
            direction=-mean_vector,
            height=float(np.cos(obj.max_angle)),
            resolution=100,
            angle=float(obj.max_angle.to("degree").magnitude),
        )
        scene.add_mesh(cone_mesh, color=cone_color, opacity=0.6)

        scene.add_scalar_bar(mapper=mapping.mapper, title="Collecting Field Real Part")

    def _add_photodiode(
        self,
        obj: Photodiode,
        scene: pv.Plotter,
        field_point_size: float = 20.0,
        cone_color: PyVistaColor = "blue") -> None:
        """
        Add a schematic photodiode disk.

        Parameters
        ----------
        obj : PyMieSim.single.detector.Photodiode
            Photodiode detector instance (unused, kept for signature uniformity).
        scene : pyvista.Plotter
            Target plotter.
        color : str or tuple, optional
            Mesh color. Defaults to "green".
        opacity : float, optional
            Mesh opacity. Defaults to 1.0.
        """
        coordinates = np.vstack(
            (
                obj._cpp_mesh.cartesian.x,
                obj._cpp_mesh.cartesian.y,
                obj._cpp_mesh.cartesian.z,
            )
        )

        points = pv.wrap(coordinates.T)

        scalar_field_real = np.asarray(obj.scalar_field).real
        abs_max = float(np.abs(scalar_field_real).max()) if scalar_field_real.size else 1.0

        mapping = scene.add_points(
            points,
            scalars=scalar_field_real,
            point_size=float(field_point_size),
            render_points_as_spheres=True,
            cmap=self.colormap,
            show_scalar_bar=False,
            clim=[-abs_max, abs_max],
        )

        mean_vector = coordinates.mean(axis=1)
        cone_mesh = pv.Cone(
            center=mean_vector / 2.0,
            direction=-mean_vector,
            height=float(np.cos(obj.max_angle)),
            resolution=100,
            angle=float(obj.max_angle.to("degree").magnitude),
        )
        scene.add_mesh(cone_mesh, color=cone_color, opacity=0.6)

        scene.add_scalar_bar(mapper=mapping.mapper, title="Collecting Field Real Part")


    def _add_data(self, data: Any, scene: pv.Plotter) -> None:
        """
        Add a `data` object that implements `_add_to_3d_ax`.

        Parameters
        ----------
        data : Any
            Object providing `_add_to_3d_ax(scene=..., colormap=...)`.
        scene : pyvista.Plotter
            Target plotter.

        Raises
        ------
        AttributeError
            If `data` does not implement `_add_to_3d_ax`.
        """
        if not hasattr(data, "_add_to_3d_ax"):
            raise AttributeError(
                "The provided `data` object cannot be added because it lacks "
                "a `_add_to_3d_ax(scene=..., colormap=...)` method."
            )

        data._add_to_3d_ax(scene=scene, colormap=self.colormap)