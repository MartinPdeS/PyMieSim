from .scatterer import Sphere, CoreShell, Cylinder  # noqa: F401
from .source import PlaneWave, Gaussian  # noqa: F401
from .detector import Photodiode, CoherentMode  # noqa: F401

import pyvista
from typing import NoReturn, Any
from MPSPlots.colormaps import blue_black_red


def plot_system(*objects: Any, data: Any = None, colormap: str = blue_black_red, show_axis_label: bool = True) -> NoReturn:
    """
    Plots a 3D visualization of a system of objects, each of which must have a `_add_to_3d_ax` method for adding itself to the scene.

    This function creates a 3D plot of multiple objects, adding each one to a shared PyVista scene.
    It ensures that each object has the required `_add_to_3d_ax` method for adding itself to the scene.
    Additionally, it includes an optional translucent sphere and the option to display axis labels.

    Args:
        objects (Any): Objects to be plotted in the 3D scene. Each object must have a `_add_to_3d_ax` method.
        show_axis_label (bool): If True, axis labels will be shown. Default is False.

    Returns:
        NoReturn: This function does not return a value. It displays the 3D visualization.
    """
    # Create a PyVista plotting scene
    scene = pyvista.Plotter(theme=pyvista.themes.DarkTheme())

    # Add each object to the scene by calling its `_add_to_3d_ax` method
    for obj in objects:
        if not hasattr(obj, '_add_to_3d_ax'):
            raise AttributeError(f'Object {obj} cannot be added to system plotting because it lacks a `_add_to_3d_ax` method.')
        obj._add_to_3d_ax(scene=scene)

    if data is not None:
        if not hasattr(data, '_add_to_3d_ax'):
            raise AttributeError(f'Data {obj} cannot be added to system plotting because it lacks a `_add_to_3d_ax` method.')
        data._add_to_3d_ax(scene=scene, colormap=colormap)
    # Add a translucent sphere to the scene to provide context or reference
    sphere = pyvista.Sphere(radius=1)
    scene.add_mesh(sphere, opacity=0.3)

    # Optionally add axis labels to the scene
    scene.add_axes_at_origin(zlabel='', xlabel='', ylabel='', labels_off=not show_axis_label)

    # Display the scene
    scene.show()
