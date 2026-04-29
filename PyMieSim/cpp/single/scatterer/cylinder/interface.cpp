#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <utility>
#include <vector>

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pint/pint.h>

#include <utils/numpy_interface.h>

#include <single/scatterer/utils.h>

#include "./cylinder.h"


namespace py = pybind11;


void register_cylinder(py::module_& module) {
    py::object ureg = get_shared_ureg();

    auto ensure_matplotlib_3d_axis = [](
        const py::object& ax
    ) -> void {
        if (!py::hasattr(ax, "get_zlim")) {
            throw std::runtime_error(
                "InfiniteCylinder plotting requires a Matplotlib 3D axis. "
                "Create it with: figure.add_subplot(111, projection=\"3d\")."
            );
        }
    };

    auto make_sphere_surface = [](
        const double radius,
        const std::size_t theta_sampling,
        const std::size_t phi_sampling
    ) -> py::tuple {
        if (theta_sampling < 2) {
            throw std::invalid_argument("theta_sampling must be at least 2.");
        }

        if (phi_sampling < 2) {
            throw std::invalid_argument("phi_sampling must be at least 2.");
        }

        constexpr double pi = 3.14159265358979323846;

        py::array_t<double> x_surface(
            py::array::ShapeContainer{
                static_cast<py::ssize_t>(phi_sampling),
                static_cast<py::ssize_t>(theta_sampling)
            }
        );

        py::array_t<double> y_surface(
            py::array::ShapeContainer{
                static_cast<py::ssize_t>(phi_sampling),
                static_cast<py::ssize_t>(theta_sampling)
            }
        );

        py::array_t<double> z_surface(
            py::array::ShapeContainer{
                static_cast<py::ssize_t>(phi_sampling),
                static_cast<py::ssize_t>(theta_sampling)
            }
        );

        auto x = x_surface.mutable_unchecked<2>();
        auto y = y_surface.mutable_unchecked<2>();
        auto z = z_surface.mutable_unchecked<2>();

        for (std::size_t phi_index = 0; phi_index < phi_sampling; ++phi_index) {
            const double phi = pi *
                static_cast<double>(phi_index) /
                static_cast<double>(phi_sampling - 1);

            for (std::size_t theta_index = 0; theta_index < theta_sampling; ++theta_index) {
                const double theta = 2.0 * pi *
                    static_cast<double>(theta_index) /
                    static_cast<double>(theta_sampling - 1);

                x(
                    static_cast<py::ssize_t>(phi_index),
                    static_cast<py::ssize_t>(theta_index)
                ) = radius * std::sin(phi) * std::cos(theta);

                y(
                    static_cast<py::ssize_t>(phi_index),
                    static_cast<py::ssize_t>(theta_index)
                ) = radius * std::sin(phi) * std::sin(theta);

                z(
                    static_cast<py::ssize_t>(phi_index),
                    static_cast<py::ssize_t>(theta_index)
                ) = radius * std::cos(phi);
            }
        }

        return py::make_tuple(x_surface, y_surface, z_surface);
    };

    auto make_cylinder_surface = [](
        const double radius,
        const double length,
        const std::size_t angular_sampling,
        const std::size_t axial_sampling
    ) -> py::tuple {
        if (angular_sampling < 3) {
            throw std::invalid_argument("angular_sampling must be at least 3.");
        }

        if (axial_sampling < 2) {
            throw std::invalid_argument("axial_sampling must be at least 2.");
        }

        constexpr double pi = 3.14159265358979323846;

        py::array_t<double> x_surface(
            py::array::ShapeContainer{
                static_cast<py::ssize_t>(axial_sampling),
                static_cast<py::ssize_t>(angular_sampling)
            }
        );

        py::array_t<double> y_surface(
            py::array::ShapeContainer{
                static_cast<py::ssize_t>(axial_sampling),
                static_cast<py::ssize_t>(angular_sampling)
            }
        );

        py::array_t<double> z_surface(
            py::array::ShapeContainer{
                static_cast<py::ssize_t>(axial_sampling),
                static_cast<py::ssize_t>(angular_sampling)
            }
        );

        auto x = x_surface.mutable_unchecked<2>();
        auto y = y_surface.mutable_unchecked<2>();
        auto z = z_surface.mutable_unchecked<2>();

        for (std::size_t axial_index = 0; axial_index < axial_sampling; ++axial_index) {
            const double y_position = -0.5 * length + length *
                static_cast<double>(axial_index) /
                static_cast<double>(axial_sampling - 1);

            for (std::size_t angular_index = 0; angular_index < angular_sampling; ++angular_index) {
                const double angle = 2.0 * pi *
                    static_cast<double>(angular_index) /
                    static_cast<double>(angular_sampling - 1);

                x(
                    static_cast<py::ssize_t>(axial_index),
                    static_cast<py::ssize_t>(angular_index)
                ) = radius * std::cos(angle);

                y(
                    static_cast<py::ssize_t>(axial_index),
                    static_cast<py::ssize_t>(angular_index)
                ) = y_position;

                z(
                    static_cast<py::ssize_t>(axial_index),
                    static_cast<py::ssize_t>(angular_index)
                ) = radius * std::sin(angle);
            }
        }

        return py::make_tuple(x_surface, y_surface, z_surface);
    };

    auto set_cylinder_axis_limits = [](
        const py::object& ax,
        const double radial_axis_limit,
        const double axial_axis_limit,
        const bool show_axes
    ) -> void {
        ax.attr("set_xlim")(-radial_axis_limit, radial_axis_limit);
        ax.attr("set_ylim")(-axial_axis_limit, axial_axis_limit);
        ax.attr("set_zlim")(-radial_axis_limit, radial_axis_limit);

        if (py::hasattr(ax, "set_box_aspect")) {
            ax.attr("set_box_aspect")(
                py::make_tuple(
                    1.0,
                    axial_axis_limit / radial_axis_limit,
                    1.0
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

    auto add_sphere_surface_to_ax = [
        make_sphere_surface
    ](
        const py::object& ax,
        const double radius,
        const py::object& color,
        const double alpha,
        const std::size_t theta_sampling,
        const std::size_t phi_sampling,
        const double linewidth,
        const bool shade
    ) -> void {
        py::tuple sphere_surface = make_sphere_surface(
            radius,
            theta_sampling,
            phi_sampling
        );

        ax.attr("plot_surface")(
            sphere_surface[0],
            sphere_surface[1],
            sphere_surface[2],
            py::arg("color") = color,
            py::arg("alpha") = alpha,
            py::arg("linewidth") = linewidth,
            py::arg("shade") = shade
        );
    };

    auto add_cylinder_surface_to_ax = [
        make_cylinder_surface
    ](
        const py::object& ax,
        const double radius,
        const double length,
        const py::object& color,
        const double alpha,
        const std::size_t angular_sampling,
        const std::size_t axial_sampling,
        const double linewidth,
        const bool shade
    ) -> void {
        py::tuple cylinder_surface = make_cylinder_surface(
            radius,
            length,
            angular_sampling,
            axial_sampling
        );

        ax.attr("plot_surface")(
            cylinder_surface[0],
            cylinder_surface[1],
            cylinder_surface[2],
            py::arg("color") = color,
            py::arg("alpha") = alpha,
            py::arg("linewidth") = linewidth,
            py::arg("shade") = shade
        );
    };

    py::class_<InfiniteCylinder, BaseScatterer, std::shared_ptr<InfiniteCylinder>>(
        module,
        "InfiniteCylinder",
        R"pbdoc(
            Infinite cylindrical scatterer for Lorenz Mie simulations.

            The cylinder is defined by its diameter, material, surrounding
            medium, and optional maximum expansion order. It exposes the
            cylindrical scattering coefficients and basic geometric properties
            used by the PyMieSim scattering engine.
        )pbdoc"
    )
        .def(
            py::init(
                [ureg](
                    const py::object& diameter,
                    const py::object& material,
                    const py::object& medium,
                    const std::size_t max_order
                ) {
                    const double diameter_meter = diameter
                        .attr("to")("meter")
                        .attr("magnitude")
                        .cast<double>();

                    const std::shared_ptr<BaseMaterial> parsed_material =
                        parse_material_object(material, ureg);

                    const std::shared_ptr<BaseMedium> parsed_medium =
                        parse_medium_object(medium, ureg);

                    return std::make_shared<InfiniteCylinder>(
                        diameter_meter,
                        std::move(parsed_material),
                        std::move(parsed_medium),
                        max_order
                    );
                }
            ),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Construct an infinite cylindrical scatterer.

                Parameters
                ----------
                diameter : pint.Quantity
                    Cylinder diameter. The value is converted internally to
                    meters.
                material : BaseMaterial or complex
                    Cylinder material model, or a constant complex refractive
                    index.
                medium : BaseMedium or float
                    Surrounding medium model, or a constant real refractive
                    index.
                max_order : int, optional
                    Maximum expansion order used in the cylindrical Mie
                    calculation. If zero, the implementation chooses an
                    appropriate truncation order.

                Notes
                -----
                The cylinder is infinite in the scattering model. The plotting
                helper draws a finite cylinder only as a normalized visual
                representation.
            )pbdoc"
        )
        .def_readonly_static(
            "property_names",
            &InfiniteCylinder::property_names,
            R"pbdoc(
                Names of the physical properties exposed by the cylinder.
            )pbdoc"
        )
        .def_property_readonly(
            "diameter",
            [ureg](
                const InfiniteCylinder& self
            ) {
                return (py::float_(self.diameter) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Cylinder diameter.

                Returns
                -------
                pint.Quantity
                    Diameter expressed as a compact length unit.
            )pbdoc"
        )
        .def(
            "print_properties",
            &InfiniteCylinder::print_properties,
            R"pbdoc(
                Print the cylinder properties to standard output.
            )pbdoc"
        )
        .def_readonly(
            "material",
            &InfiniteCylinder::material,
            R"pbdoc(
                Material model of the cylinder.
            )pbdoc"
        )
        .def_readonly(
            "medium",
            &InfiniteCylinder::medium,
            R"pbdoc(
                Medium surrounding the cylinder.
            )pbdoc"
        )
        .def_property_readonly(
            "radius",
            [ureg](
                const InfiniteCylinder& self
            ) {
                return (py::float_(self.diameter / 2.0) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Cylinder radius.

                Returns
                -------
                pint.Quantity
                    Radius expressed as a compact length unit.
            )pbdoc"
        )
        .def_property_readonly(
            "a1n",
            [](
                InfiniteCylinder& self
            ) {
                return vector_as_numpy_view(self, self.a1n);
            },
            R"pbdoc(
                Cylindrical scattering coefficients \(a_{1n}\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the \(a_{1n}\)
                    coefficients.
            )pbdoc"
        )
        .def_property_readonly(
            "b1n",
            [](
                InfiniteCylinder& self
            ) {
                return vector_as_numpy_view(self, self.b1n);
            },
            R"pbdoc(
                Cylindrical scattering coefficients \(b_{1n}\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the \(b_{1n}\)
                    coefficients.
            )pbdoc"
        )
        .def_property_readonly(
            "a2n",
            [](
                InfiniteCylinder& self
            ) {
                return vector_as_numpy_view(self, self.a2n);
            },
            R"pbdoc(
                Cylindrical scattering coefficients \(a_{2n}\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the \(a_{2n}\)
                    coefficients.
            )pbdoc"
        )
        .def_property_readonly(
            "b2n",
            [](
                InfiniteCylinder& self
            ) {
                return vector_as_numpy_view(self, self.b2n);
            },
            R"pbdoc(
                Cylindrical scattering coefficients \(b_{2n}\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the \(b_{2n}\)
                    coefficients.
            )pbdoc"
        )
        .def(
            "_add_to_ax",
            [
                ensure_matplotlib_3d_axis,
                add_cylinder_surface_to_ax,
                add_sphere_surface_to_ax,
                set_cylinder_axis_limits
            ](
                const InfiniteCylinder& self,
                const py::object& ax,
                const py::object& color,
                const double alpha,
                const double display_radius,
                const double display_length,
                const bool show_unit_sphere,
                const py::object& unit_sphere_color,
                const double unit_sphere_alpha,
                const std::size_t angular_sampling,
                const std::size_t axial_sampling,
                const std::size_t unit_sphere_theta_sampling,
                const std::size_t unit_sphere_phi_sampling,
                const double linewidth,
                const bool shade,
                const bool show_axes
            ) -> void {
                ensure_matplotlib_3d_axis(ax);

                add_cylinder_surface_to_ax(
                    ax,
                    display_radius,
                    display_length,
                    color,
                    alpha,
                    angular_sampling,
                    axial_sampling,
                    linewidth,
                    shade
                );

                if (show_unit_sphere) {
                    add_sphere_surface_to_ax(
                        ax,
                        1.0,
                        unit_sphere_color,
                        unit_sphere_alpha,
                        unit_sphere_theta_sampling,
                        unit_sphere_phi_sampling,
                        linewidth,
                        shade
                    );
                }

                const double radial_axis_limit = std::max(1.0, display_radius);
                const double axial_axis_limit = std::max(1.0, display_length / 2.0);

                set_cylinder_axis_limits(
                    ax,
                    radial_axis_limit,
                    axial_axis_limit,
                    show_axes
                );
            },
            py::arg("ax"),
            py::arg("color") = py::str("black"),
            py::arg("alpha") = 1.0,
            py::arg("display_radius") = 0.1,
            py::arg("display_length") = 2.0,
            py::arg("show_unit_sphere") = true,
            py::arg("unit_sphere_color") = py::str("grey"),
            py::arg("unit_sphere_alpha") = 0.15,
            py::arg("angular_sampling") = 48,
            py::arg("axial_sampling") = 8,
            py::arg("unit_sphere_theta_sampling") = 24,
            py::arg("unit_sphere_phi_sampling") = 24,
            py::arg("linewidth") = 0.0,
            py::arg("shade") = false,
            py::arg("show_axes") = false,
            R"pbdoc(
                Add the infinite cylinder geometry to a Matplotlib 3D axis.

                This private plotting helper replaces the previous PyVista based
                ``add_to_scene`` method. It draws a finite normalized cylinder
                as a visual representation of the infinite cylindrical scatterer
                and, optionally, a transparent unit sphere for angular reference.

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Matplotlib axis created with ``projection="3d"``.
                color : str or tuple, optional
                    Matplotlib color used for the displayed cylinder.
                alpha : float, optional
                    Transparency of the displayed cylinder.
                display_radius : float, optional
                    Radius used for visualization. This is a display scale, not
                    the physical radius of the particle.
                display_length : float, optional
                    Finite displayed cylinder length. This is only a
                    visualization length because the physical model is an
                    infinite cylinder.
                show_unit_sphere : bool, optional
                    If ``True``, draw a transparent unit sphere as an angular
                    reference surface.
                unit_sphere_color : str or tuple, optional
                    Matplotlib color used for the unit sphere.
                unit_sphere_alpha : float, optional
                    Transparency of the unit sphere.
                angular_sampling : int, optional
                    Number of angular samples used to draw the cylinder surface.
                axial_sampling : int, optional
                    Number of axial samples used to draw the cylinder surface.
                unit_sphere_theta_sampling : int, optional
                    Number of azimuthal samples used to draw the unit sphere.
                unit_sphere_phi_sampling : int, optional
                    Number of polar samples used to draw the unit sphere.
                linewidth : float, optional
                    Surface mesh line width passed to Matplotlib.
                shade : bool, optional
                    If ``True``, enable Matplotlib surface shading.
                show_axes : bool, optional
                    If ``True``, display Cartesian axis labels, ticks, and panes.
                    If ``False``, hide the Matplotlib 3D axis frame after setting
                    the plotting limits.

                Notes
                -----
                The plotted cylinder has a finite length only for visualization.
                The scattering calculation still represents an infinite cylinder.
                The cylinder is drawn along the y axis to match the previous
                PyVista based visualization.
            )pbdoc"
        )
    ;
}