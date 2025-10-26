#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import MPSPlots
from pydantic.dataclasses import dataclass
import matplotlib.pyplot as plt
from TypedUnit import ureg

from PyMieSim.utils import config_dict, rotate_on_x


@dataclass(config=config_dict, kw_only=True)
class Footprint:
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
            -max_angle.to("radian")
            .magnitude : max_angle.to("radian")
            .magnitude : n_point,
            0 : numpy.pi : n_point,
        ]

        max_distance_direct_space = 1 / (
            numpy.sin(max_angle) * self.scatterer.source.wavenumber / (2 * numpy.pi)
        )

        x = y = (
            numpy.linspace(-1, 1, self.sampling)
            * self.sampling
            / 2
            * max_distance_direct_space
            / self.padding_factor
        )

        _, phi, theta = rotate_on_x(phi + numpy.pi / 2, theta, numpy.pi / 2)

        far_field_para, far_field_perp = self.scatterer.get_farfields_array(
            phi=(phi.ravel() + numpy.pi / 2) * ureg.radian,
            theta=theta.ravel() * ureg.radian,
            r=1.0 * ureg.meter,
        )

        detector_structured_farfield = self.detector.get_structured_scalarfield(
            sampling=self.sampling
        )

        perpendicular_projection = (
            detector_structured_farfield * far_field_perp.reshape(theta.shape)
        )

        parallel_projection = detector_structured_farfield * far_field_para.reshape(
            theta.shape
        )

        fourier_parallel = self.get_fourier_component(parallel_projection)
        fourier_perpendicular = self.get_fourier_component(perpendicular_projection)

        self.mapping = fourier_parallel + fourier_perpendicular
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
        padded_scalar = numpy.pad(
            scalar,
            pad_width=((offset, offset), (offset, offset)),
            mode="constant",
            constant_values=0,
        )

        # Perform the two-dimensional inverse Fourier transform on the padded scalar field.
        fourier_transformed = numpy.fft.ifft2(padded_scalar)

        # Compute the squared magnitude and center the zero-frequency component.
        fourier_magnitude_squared = (
            numpy.abs(numpy.fft.fftshift(fourier_transformed)) ** 2
        )

        # Extract the central portion corresponding to the original sampling rate adjusted by the padding factor.
        central_portion = fourier_magnitude_squared[offset:-offset, offset:-offset]

        return central_portion

    def plot(self, colormap: str = "gray") -> None:
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
                title="Scatterer Footprint",
                xlabel=r"Offset distance in X-axis [$\mu$m]",
                ylabel=r"Offset distance in Y-axis [$\mu$m]",
            )

            ax.pcolormesh(
                self.direct_y.magnitude,
                self.direct_x.magnitude,
                self.mapping,
                cmap=colormap,
            )

            plt.show()
