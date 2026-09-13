"""Common array conversion and rendering helpers for angular representations."""

import numpy
import matplotlib
from matplotlib import colors
from matplotlib.colors import Colormap


def resolve_colormap(colormap: str | Colormap | None) -> Colormap:
    """Resolve a registered name or return the supplied colormap."""
    return matplotlib.colormaps.get_cmap(colormap)


def signed_normalization(field: numpy.ndarray, percentile_clip: float | None) -> colors.Normalize:
    """Use symmetric limits, including for zero or non-finite fields."""
    finite = numpy.abs(field[numpy.isfinite(field)])
    maximum = 1.0 if finite.size == 0 else (
        numpy.max(finite) if percentile_clip is None else numpy.percentile(finite, percentile_clip)
    )
    if not numpy.isfinite(maximum) or maximum <= 0:
        maximum = 1.0
    return colors.Normalize(vmin=-maximum, vmax=maximum)


def intensity_normalization(
    field: numpy.ndarray,
    scale: str,
    percentile_clip: float | None,
    *,
    include_zeros: bool = False,
) -> colors.Normalize:
    """Normalize intensity, retaining the SPF convention of including zeros in clipping."""
    finite = field[numpy.isfinite(field)]
    values = finite[finite >= 0] if include_zeros else finite[finite > 0]
    if values.size == 0:
        return colors.Normalize(vmin=0, vmax=1)
    upper = numpy.max(values) if percentile_clip is None else numpy.percentile(values, percentile_clip)
    if not numpy.isfinite(upper) or upper <= 0:
        upper = 1.0
    if scale == "linear":
        return colors.Normalize(vmin=0, vmax=upper)
    if scale == "log":
        positive = values[values > 0]
        if positive.size == 0:
            return colors.Normalize(vmin=0, vmax=upper)
        lower = max(numpy.max(positive) * 1e-6, numpy.min(positive))
        if lower >= upper:
            lower = upper * 1e-6
        return colors.LogNorm(vmin=lower, vmax=upper)
    raise ValueError("Invalid intensity scale. Expected 'linear' or 'log'.")


class AngularPlotMixin:
    """Rendering operations shared by representations sampled on a square mesh."""

    sampling: int

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
                "Cannot reshape array to the structured angular mesh. "
                f"Expected {expected_size} values from sampling={self.sampling}, "
                f"but received {flat_array.size}."
            )

        return flat_array.reshape(
            (self.sampling, self.sampling),
            order="F",
        )
