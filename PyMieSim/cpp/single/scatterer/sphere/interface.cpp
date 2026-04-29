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

#include "./sphere.h"


namespace py = pybind11;


void register_sphere(py::module_& module) {
    py::object ureg = get_shared_ureg();

    auto ensure_matplotlib_3d_axis = [](
        const py::object& ax
    ) -> void {
        if (!py::hasattr(ax, "get_zlim")) {
            throw std::runtime_error(
                "Sphere plotting requires a Matplotlib 3D axis. "
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

    auto set_sphere_axis_limits = [](
        const py::object& ax,
        const double axis_limit,
        const bool show_axes
    ) -> void {
        ax.attr("set_xlim")(-axis_limit, axis_limit);
        ax.attr("set_ylim")(-axis_limit, axis_limit);
        ax.attr("set_zlim")(-axis_limit, axis_limit);

        if (py::hasattr(ax, "set_box_aspect")) {
            ax.attr("set_box_aspect")(py::make_tuple(1.0, 1.0, 1.0));
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

    py::class_<Sphere, BaseScatterer, std::shared_ptr<Sphere>>(
        module,
        "Sphere",
        R"pbdoc(
            Spherical scatterer for Lorenz Mie simulations.

            The sphere is defined by its diameter, material, surrounding medium,
            and optional maximum multipole order. It exposes the Mie scattering
            coefficients and basic geometric properties used by the PyMieSim
            scattering engine.
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
                        .attr("to")(ureg.attr("meter"))
                        .attr("magnitude")
                        .cast<double>();

                    const std::shared_ptr<BaseMaterial> parsed_material =
                        parse_material_object(material, ureg);

                    const std::shared_ptr<BaseMedium> parsed_medium =
                        parse_medium_object(medium, ureg);

                    return std::make_shared<Sphere>(
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
                Construct a spherical scatterer.

                Parameters
                ----------
                diameter : pint.Quantity
                    Sphere diameter. The value is converted internally to meters.
                material : BaseMaterial or complex
                    Sphere material model, or a constant complex refractive index.
                medium : BaseMedium or float
                    Surrounding medium model, or a constant real refractive index.
                max_order : int, optional
                    Maximum multipole order used in the Mie expansion. If zero,
                    the implementation chooses an appropriate truncation order.

                Notes
                -----
                The material and medium may be dispersive objects or constant
                refractive indices, depending on the Python object passed to the
                constructor.
            )pbdoc"
        )
        .def_readonly_static(
            "property_names",
            &Sphere::property_names,
            R"pbdoc(
                Names of the physical properties exposed by the sphere.
            )pbdoc"
        )
        .def(
            "print_properties",
            &Sphere::print_properties,
            R"pbdoc(
                Print the sphere properties to standard output.
            )pbdoc"
        )
        .def_property_readonly(
            "diameter",
            [ureg](
                const Sphere& self
            ) {
                return (py::float_(self.diameter) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Sphere diameter.

                Returns
                -------
                pint.Quantity
                    Diameter expressed as a compact length unit.
            )pbdoc"
        )
        .def_readonly(
            "material",
            &Sphere::material,
            R"pbdoc(
                Material model of the sphere.
            )pbdoc"
        )
        .def_readonly(
            "medium",
            &Sphere::medium,
            R"pbdoc(
                Medium surrounding the sphere.
            )pbdoc"
        )
        .def_property_readonly(
            "radius",
            [ureg](
                const Sphere& self
            ) {
                return (py::float_(self.diameter / 2.0) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Sphere radius.

                Returns
                -------
                pint.Quantity
                    Radius expressed as a compact length unit.
            )pbdoc"
        )
        .def_property_readonly(
            "volume",
            [ureg](
                const Sphere& self
            ) {
                const double radius = self.diameter / 2.0;
                const double volume = (4.0 / 3.0) * Constants::PI * std::pow(radius, 3);

                return (py::float_(volume) * ureg.attr("meter**3")).attr("to_compact")();
            },
            R"pbdoc(
                Sphere volume.

                Returns
                -------
                pint.Quantity
                    Volume expressed as a compact cubic length unit.
            )pbdoc"
        )
        .def_property_readonly(
            "an",
            [](
                Sphere& self
            ) {
                return vector_as_numpy_view(self, self.an);
            },
            R"pbdoc(
                Electric type external Mie coefficients \(a_n\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the external scattering
                    coefficients \(a_n\).
            )pbdoc"
        )
        .def_property_readonly(
            "bn",
            [](
                Sphere& self
            ) {
                return vector_as_numpy_view(self, self.bn);
            },
            R"pbdoc(
                Magnetic type external Mie coefficients \(b_n\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the external scattering
                    coefficients \(b_n\).
            )pbdoc"
        )
        .def_property_readonly(
            "cn",
            [](
                Sphere& self
            ) {
                return vector_as_numpy_view(self, self.cn);
            },
            R"pbdoc(
                Electric type internal Mie coefficients \(c_n\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the internal coefficients
                    \(c_n\).
            )pbdoc"
        )
        .def_property_readonly(
            "dn",
            [](
                Sphere& self
            ) {
                return vector_as_numpy_view(self, self.dn);
            },
            R"pbdoc(
                Magnetic type internal Mie coefficients \(d_n\).

                Returns
                -------
                numpy.ndarray
                    Zero copy complex array containing the internal coefficients
                    \(d_n\).
            )pbdoc"
        )
        .def(
            "_add_to_ax",
            [
                ensure_matplotlib_3d_axis,
                add_sphere_surface_to_ax,
                set_sphere_axis_limits
            ](
                const Sphere& self,
                const py::object& ax,
                const py::object& color,
                const double alpha,
                const double display_radius,
                const bool show_unit_sphere,
                const py::object& unit_sphere_color,
                const double unit_sphere_alpha,
                const std::size_t theta_sampling,
                const std::size_t phi_sampling,
                const double linewidth,
                const bool shade,
                const bool show_axes
            ) -> void {
                ensure_matplotlib_3d_axis(ax);

                add_sphere_surface_to_ax(
                    ax,
                    display_radius,
                    color,
                    alpha,
                    theta_sampling,
                    phi_sampling,
                    linewidth,
                    shade
                );

                if (show_unit_sphere) {
                    add_sphere_surface_to_ax(
                        ax,
                        1.0,
                        unit_sphere_color,
                        unit_sphere_alpha,
                        theta_sampling,
                        phi_sampling,
                        linewidth,
                        shade
                    );
                }

                const double axis_limit = std::max(1.0, display_radius);

                set_sphere_axis_limits(
                    ax,
                    axis_limit,
                    show_axes
                );
            },
            py::arg("ax"),
            py::arg("color") = py::str("black"),
            py::arg("alpha") = 1.0,
            py::arg("display_radius") = 0.1,
            py::arg("show_unit_sphere") = true,
            py::arg("unit_sphere_color") = py::str("grey"),
            py::arg("unit_sphere_alpha") = 0.15,
            py::arg("theta_sampling") = 24,
            py::arg("phi_sampling") = 24,
            py::arg("linewidth") = 0.0,
            py::arg("shade") = false,
            py::arg("show_axes") = false,
            R"pbdoc(
                Add the sphere geometry to a Matplotlib 3D axis.

                This private plotting helper replaces the previous PyVista based
                ``add_to_scene`` method. It draws a normalized visual sphere and,
                optionally, a transparent unit sphere for angular reference.

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Matplotlib axis created with ``projection="3d"``.
                color : str or tuple, optional
                    Matplotlib color used for the displayed sphere.
                alpha : float, optional
                    Transparency of the displayed sphere.
                display_radius : float, optional
                    Radius used for visualization. This is a display scale, not
                    the physical radius of the particle.
                show_unit_sphere : bool, optional
                    If ``True``, draw a transparent unit sphere as an angular
                    reference surface.
                unit_sphere_color : str or tuple, optional
                    Matplotlib color used for the unit sphere.
                unit_sphere_alpha : float, optional
                    Transparency of the unit sphere.
                theta_sampling : int, optional
                    Number of azimuthal samples used to draw the surfaces.
                phi_sampling : int, optional
                    Number of polar samples used to draw the surfaces.
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
                The physical particle diameter is usually many orders of
                magnitude smaller than the angular reference sphere used in
                scattering visualizations. For that reason, ``display_radius``
                controls only the plotted size.
            )pbdoc"
        )
        ;
}