#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <pint/pint.h>

#include <utils/casting.h>
#include <utils/numpy_interface.h>

#include "./source.h"


namespace py = pybind11;


PYBIND11_MODULE(source, module)
{
    py::object ureg = get_shared_ureg();

    module.doc() = R"pbdoc(
        Source bindings for PyMieSim.

        This module exposes optical excitation sources used by the single
        scatterer simulation API. Plotting helpers use Matplotlib 3D axes through
        private ``_add_to_ax`` methods. Figure creation, styling, and public
        visualization APIs should remain controlled by the Python layer.
    )pbdoc";

    auto ensure_matplotlib_3d_axis = [](
        const py::object& ax
    ) -> void {
        if (!py::hasattr(ax, "get_zlim")) {
            throw std::runtime_error(
                "Source plotting requires a Matplotlib 3D axis. "
                "Create it with: figure.add_subplot(111, projection=\"3d\")."
            );
        }
    };

    auto set_source_axis_limits = [](
        const py::object& ax,
        const double radial_axis_limit,
        const double axial_axis_limit,
        const bool show_axes
    ) -> void {
        ax.attr("set_xlim")(-radial_axis_limit, radial_axis_limit);
        ax.attr("set_ylim")(-radial_axis_limit, radial_axis_limit);
        ax.attr("set_zlim")(-axial_axis_limit, axial_axis_limit);

        if (py::hasattr(ax, "set_box_aspect")) {
            ax.attr("set_box_aspect")(
                py::make_tuple(
                    1.0,
                    1.0,
                    axial_axis_limit / radial_axis_limit
                )
            );
        }

        if (show_axes) {
            ax.attr("set_xlabel")("x");
            ax.attr("set_ylabel")("y");
            ax.attr("set_zlabel")("z");
        } else {
            ax.attr("set_axis_off")();
        }
    };

    auto add_vector_arrow_to_ax = [](
        const py::object& ax,
        const std::array<double, 3>& origin,
        const std::array<double, 3>& direction,
        const py::object& color,
        const double length,
        const double arrow_length_ratio,
        const double linewidth,
        const double alpha
    ) -> void {
        ax.attr("quiver")(
            py::float_(origin[0]),
            py::float_(origin[1]),
            py::float_(origin[2]),
            py::float_(direction[0]),
            py::float_(direction[1]),
            py::float_(direction[2]),
            py::arg("length") = length,
            py::arg("arrow_length_ratio") = arrow_length_ratio,
            py::arg("linewidth") = linewidth,
            py::arg("color") = color,
            py::arg("alpha") = alpha,
            py::arg("normalize") = true
        );
    };

    auto add_source_vectors_to_ax = [
        add_vector_arrow_to_ax
    ](
        const BaseSource& source,
        const py::object& ax,
        const py::object& propagation_color,
        const py::object& electric_field_color,
        const double origin_z,
        const double propagation_length,
        const double electric_field_length,
        const double arrow_length_ratio,
        const double linewidth,
        const double alpha
    ) -> void {
        const std::complex<double> electric_field_x = source.polarization.jones_vector[0];
        const std::complex<double> electric_field_y = source.polarization.jones_vector[1];

        const double electric_x = std::real(electric_field_x);
        const double electric_y = std::real(electric_field_y);

        const double electric_norm = std::sqrt(
            electric_x * electric_x +
            electric_y * electric_y
        );

        std::array<double, 3> normalized_electric_direction = {1.0, 0.0, 0.0};

        if (electric_norm > 0.0) {
            normalized_electric_direction = {
                electric_x / electric_norm,
                electric_y / electric_norm,
                0.0
            };
        }

        const std::array<double, 3> origin = {0.0, 0.0, origin_z};
        const std::array<double, 3> propagation_direction = {0.0, 0.0, 1.0};

        add_vector_arrow_to_ax(
            ax,
            origin,
            propagation_direction,
            propagation_color,
            propagation_length,
            arrow_length_ratio,
            linewidth,
            alpha
        );

        add_vector_arrow_to_ax(
            ax,
            origin,
            normalized_electric_direction,
            electric_field_color,
            electric_field_length,
            arrow_length_ratio,
            linewidth,
            alpha
        );
    };

    py::class_<BaseSource, std::shared_ptr<BaseSource>>(
        module,
        "BaseSource",
        R"pbdoc(
            Base class for optical excitation sources.

            A ``BaseSource`` defines the incident field that illuminates a
            scatterer. Concrete sources provide the wavelength, electric field
            amplitude, and polarization state.

            Notes
            -----
            Units are returned as Pint quantities using the shared unit registry.
            Polarization is represented internally as a Jones vector with complex
            values.
        )pbdoc"
    )
        .def(
            py::init<>(),
            R"pbdoc(
                Construct a source base instance.

                Notes
                -----
                This constructor is primarily exposed for binding completeness.
                In user code, instantiate a concrete source such as
                :class:`~PyMieSim.single.source.PlaneWave` or
                :class:`~PyMieSim.single.source.Gaussian`.
            )pbdoc"
        )
        .def_property_readonly(
            "wavelength",
            [ureg](
                const BaseSource& self
            ) {
                return py::float_(self.wavelength) * ureg.attr("meter");
            },
            R"pbdoc(
                Wavelength in vacuum.

                Returns
                -------
                pint.Quantity
                    Wavelength with units of length.
            )pbdoc"
        )
        .def_property_readonly(
            "wavenumber_vacuum",
            [ureg](
                const BaseSource& self
            ) {
                return (py::float_(self.wavenumber_vacuum) * ureg.attr("1/meter")).attr("to_compact")();
            },
            R"pbdoc(
                Vacuum wavenumber.

                This is defined as ``k0 = 2*pi / wavelength``.

                Returns
                -------
                pint.Quantity
                    Wavenumber with units of inverse length.
            )pbdoc"
        )
        .def_property_readonly(
            "amplitude",
            [ureg](
                const BaseSource& self
            ) {
                return (py::float_(self.amplitude) * ureg.attr("volt/meter")).attr("to_compact")();
            },
            R"pbdoc(
                Electric field amplitude.

                Returns
                -------
                pint.Quantity
                    Field amplitude with units of electric field.
            )pbdoc"
        )
        .def_readonly(
            "polarization",
            &BaseSource::polarization,
            R"pbdoc(
                Polarization state of the source.

                Returns
                -------
                PolarizationState
                    Polarization state represented internally by a Jones vector.
            )pbdoc"
        )
        .def(
            "_add_to_ax",
            [
                ensure_matplotlib_3d_axis,
                add_source_vectors_to_ax,
                set_source_axis_limits
            ](
                const BaseSource& self,
                const py::object& ax,
                const py::object& propagation_color,
                const py::object& electric_field_color,
                const double origin_z,
                const double propagation_length,
                const double electric_field_length,
                const double arrow_length_ratio,
                const double linewidth,
                const double alpha,
                const bool show_axes
            ) -> void {
                ensure_matplotlib_3d_axis(ax);

                add_source_vectors_to_ax(
                    self,
                    ax,
                    propagation_color,
                    electric_field_color,
                    origin_z,
                    propagation_length,
                    electric_field_length,
                    arrow_length_ratio,
                    linewidth,
                    alpha
                );

                const double radial_axis_limit = std::max(1.0, electric_field_length);
                const double axial_axis_limit = std::max(
                    1.0,
                    std::abs(origin_z) + propagation_length
                );

                set_source_axis_limits(
                    ax,
                    radial_axis_limit,
                    axial_axis_limit,
                    show_axes
                );
            },
            py::arg("ax"),
            py::arg("propagation_color") = py::str("red"),
            py::arg("electric_field_color") = py::str("blue"),
            py::arg("origin_z") = -2.1,
            py::arg("propagation_length") = 1.1,
            py::arg("electric_field_length") = 0.8,
            py::arg("arrow_length_ratio") = 0.18,
            py::arg("linewidth") = 4.0,
            py::arg("alpha") = 1.0,
            py::arg("show_axes") = false,
            R"pbdoc(
                Add the source vectors to a Matplotlib 3D axis.

                This private plotting helper replaces the previous PyVista based
                ``add_to_scene`` method. It draws the propagation vector and the
                real part of the electric field polarization vector.

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Matplotlib axis created with ``projection="3d"``.
                propagation_color : str or tuple, optional
                    Matplotlib color used for the propagation direction.
                electric_field_color : str or tuple, optional
                    Matplotlib color used for the electric field vector.
                origin_z : float, optional
                    Starting z coordinate for both arrows.
                propagation_length : float, optional
                    Display length of the propagation vector.
                electric_field_length : float, optional
                    Display length of the electric field vector.
                arrow_length_ratio : float, optional
                    Fraction of the arrow length occupied by the arrow head.
                linewidth : float, optional
                    Arrow line width.
                alpha : float, optional
                    Arrow transparency.
                show_axes : bool, optional
                    If ``True``, display Cartesian axis labels, ticks, and panes.
                    If ``False``, hide the Matplotlib 3D axis frame after setting
                    the plotting limits.

                Notes
                -----
                This method is intended for visualization only and does not affect
                the physical properties of the source.
            )pbdoc"
        )
        ;

    py::class_<Planewave, BaseSource, std::shared_ptr<Planewave>>(
        module,
        "PlaneWave",
        R"pbdoc(
            Monochromatic plane wave source.

            This source models a spatially uniform incident field with a
            specified wavelength, polarization, and electric field amplitude.
        )pbdoc"
    )
        .def(
            py::init(
                [](
                    const py::object& wavelength,
                    const py::object& polarization,
                    const py::object& amplitude
                ) {
                    return std::make_shared<Planewave>(
                        wavelength.attr("to")("meter").attr("magnitude").cast<double>(),
                        Casting::Polarization::cast_py_to_polarization_state(polarization),
                        amplitude.attr("to")("volt/meter").attr("magnitude").cast<double>()
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("amplitude"),
            R"pbdoc(
                Construct a plane wave source.

                Parameters
                ----------
                wavelength : pint.Quantity
                    Wavelength in vacuum. Must be convertible to meters.
                polarization : object
                    Polarization specification. Accepted inputs include bound
                    C++ Jones vector objects, PyMieSim polarization objects, or
                    any object accepted by the linear polarization constructor.
                amplitude : pint.Quantity
                    Electric field amplitude. Must be convertible to V/m.
            )pbdoc"
        )
        ;

    py::class_<Gaussian, BaseSource, std::shared_ptr<Gaussian>>(
        module,
        "Gaussian",
        R"pbdoc(
            Monochromatic Gaussian beam source.

            This source models a diffraction limited Gaussian beam characterized
            by wavelength, polarization, numerical aperture, and optical power.
        )pbdoc"
    )
        .def(
            py::init(
                [](
                    const py::object& wavelength,
                    const py::object& polarization,
                    const py::object& numerical_aperture,
                    const py::object& optical_power
                ) {
                    return std::make_shared<Gaussian>(
                        wavelength.attr("to")("meter").attr("magnitude").cast<double>(),
                        Casting::Polarization::cast_py_to_polarization_state(polarization),
                        numerical_aperture.cast<double>(),
                        optical_power.attr("to")("watt").attr("magnitude").cast<double>()
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("numerical_aperture"),
            py::arg("optical_power"),
            R"pbdoc(
                Construct a Gaussian beam source.

                Parameters
                ----------
                wavelength : pint.Quantity
                    Wavelength in vacuum. Must be convertible to meters.
                polarization : object
                    Polarization specification. See
                    :class:`~PyMieSim.single.source.PlaneWave` for accepted
                    formats.
                numerical_aperture : float
                    Dimensionless numerical aperture.
                optical_power : pint.Quantity
                    Optical power. Must be convertible to watts.
            )pbdoc"
        )
        .def_property_readonly(
            "waist",
            [ureg](
                const Gaussian& self
            ) {
                return (py::float_(self.waist) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Beam waist radius.

                Returns
                -------
                pint.Quantity
                    Waist radius with units of length.
            )pbdoc"
        )
        .def_property_readonly(
            "peak_intensity",
            [ureg](
                const Gaussian& self
            ) {
                return (py::float_(self.peak_intensity) * ureg.attr("watt/meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Peak on axis intensity.

                Returns
                -------
                pint.Quantity
                    Peak intensity with units of power per area.
            )pbdoc"
        )
        .def_property_readonly(
            "area",
            [ureg](
                const Gaussian& self
            ) {
                return (py::float_(self.area) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Effective beam area.

                Returns
                -------
                pint.Quantity
                    Area with units of square length.
            )pbdoc"
        )
        .def_property_readonly(
            "numerical_aperture",
            [](
                const Gaussian& self
            ) {
                return py::float_(self.numerical_aperture);
            },
            R"pbdoc(
                Numerical aperture.

                Returns
                -------
                float
                    Dimensionless numerical aperture.
            )pbdoc"
        )
        .def_property_readonly(
            "optical_power",
            [ureg](
                const Gaussian& self
            ) {
                return (py::float_(self.optical_power) * ureg.attr("watt")).attr("to_compact")();
            },
            R"pbdoc(
                Optical power.

                Returns
                -------
                pint.Quantity
                    Optical power with units of watts.
            )pbdoc"
        )
        ;
}