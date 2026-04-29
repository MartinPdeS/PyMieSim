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

#include "./coreshell.h"


namespace py = pybind11;


void register_coreshell(py::module_& module) {
    py::object ureg = get_shared_ureg();

    auto ensure_matplotlib_3d_axis = [](
        const py::object& ax
    ) -> void {
        if (!py::hasattr(ax, "get_zlim")) {
            throw std::runtime_error(
                "CoreShell plotting requires a Matplotlib 3D axis. "
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

    auto set_coreshell_axis_limits = [](
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

    py::class_<CoreShell, BaseScatterer, std::shared_ptr<CoreShell>>(
        module,
        "CoreShell",
        R"pbdoc(
            Core shell spherical scatterer for Lorenz Mie simulations.

            The particle is defined by a spherical core, a concentric shell,
            their respective optical materials, and the surrounding medium. It
            exposes the Mie coefficients and geometric properties required by the
            PyMieSim scattering engine.
        )pbdoc"
    )
        .def(
            py::init(
                [ureg](
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_material,
                    const py::object& shell_material,
                    const py::object& medium,
                    const int max_order
                ) {
                    const double core_diameter_meter = core_diameter
                        .attr("to")(ureg.attr("meter"))
                        .attr("magnitude")
                        .cast<double>();

                    const double shell_thickness_meter = shell_thickness
                        .attr("to")(ureg.attr("meter"))
                        .attr("magnitude")
                        .cast<double>();

                    const std::shared_ptr<BaseMaterial> parsed_core_material =
                        parse_material_object(core_material, ureg);

                    const std::shared_ptr<BaseMaterial> parsed_shell_material =
                        parse_material_object(shell_material, ureg);

                    const std::shared_ptr<BaseMedium> parsed_medium =
                        parse_medium_object(medium, ureg);

                    return std::make_shared<CoreShell>(
                        core_diameter_meter,
                        shell_thickness_meter,
                        std::move(parsed_core_material),
                        std::move(parsed_shell_material),
                        std::move(parsed_medium),
                        max_order
                    );
                }
            ),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_material"),
            py::arg("shell_material"),
            py::arg("medium"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Construct a core shell spherical scatterer.

                Parameters
                ----------
                core_diameter : pint.Quantity
                    Diameter of the spherical core. The value is converted
                    internally to meters.
                shell_thickness : pint.Quantity
                    Thickness of the shell surrounding the core. The value is
                    converted internally to meters.
                core_material : BaseMaterial or complex
                    Core material model, or a constant complex refractive index.
                shell_material : BaseMaterial or complex
                    Shell material model, or a constant complex refractive index.
                medium : BaseMedium or float
                    Surrounding medium model, or a constant real refractive index.
                max_order : int, optional
                    Maximum multipole order used in the Mie expansion. If zero,
                    the implementation chooses an appropriate truncation order.

                Notes
                -----
                The external particle radius is
                ``core_diameter / 2 + shell_thickness``. The plotting helper uses
                normalized display radii and does not draw the particle at its
                physical nanometer scale.
            )pbdoc"
        )
        .def(
            "print_properties",
            &CoreShell::print_properties,
            R"pbdoc(
                Print the core shell scatterer properties to standard output.
            )pbdoc"
        )
        .def_property_readonly(
            "core_diameter",
            [ureg](
                const CoreShell& self
            ) {
                return (py::float_(self.core_diameter) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Diameter of the core.

                Returns
                -------
                pint.Quantity
                    Core diameter expressed as a compact length unit.
            )pbdoc"
        )
        .def_property_readonly(
            "shell_thickness",
            [ureg](
                const CoreShell& self
            ) {
                return (py::float_(self.shell_thickness) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Thickness of the shell.

                Returns
                -------
                pint.Quantity
                    Shell thickness expressed as a compact length unit.
            )pbdoc"
        )
        .def_readonly(
            "core_material",
            &CoreShell::core_material,
            R"pbdoc(
                Material model of the core.
            )pbdoc"
        )
        .def_readonly(
            "shell_material",
            &CoreShell::shell_material,
            R"pbdoc(
                Material model of the shell.
            )pbdoc"
        )
        .def_readonly_static(
            "property_names",
            &CoreShell::property_names,
            R"pbdoc(
                Names of the physical properties exposed by the core shell scatterer.
            )pbdoc"
        )
        .def_property_readonly(
            "radius",
            [ureg](
                const CoreShell& self
            ) {
                const double radius = self.core_diameter / 2.0 + self.shell_thickness;

                return (py::float_(radius) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Outer radius of the core shell particle.

                Returns
                -------
                pint.Quantity
                    Outer radius expressed as a compact length unit.
            )pbdoc"
        )
        .def_property_readonly(
            "volume",
            [ureg](
                const CoreShell& self
            ) {
                const double core_radius = self.core_diameter / 2.0;
                const double shell_outer_radius = core_radius + self.shell_thickness;

                const double volume = (4.0 / 3.0) * Constants::PI * (
                    shell_outer_radius * shell_outer_radius * shell_outer_radius -
                    core_radius * core_radius * core_radius
                );

                return (py::float_(volume) * ureg.attr("meter**3")).attr("to_compact")();
            },
            R"pbdoc(
                Shell volume excluding the core volume.

                Returns
                -------
                pint.Quantity
                    Shell volume expressed as a compact cubic length unit.
            )pbdoc"
        )
        .def_property_readonly(
            "an",
            [](
                CoreShell& self
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
                CoreShell& self
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
                CoreShell& self
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
                CoreShell& self
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
                set_coreshell_axis_limits
            ](
                const CoreShell& self,
                const py::object& ax,
                const py::object& core_color,
                const py::object& shell_color,
                const double core_alpha,
                const double shell_alpha,
                const double display_radius,
                const bool show_core,
                const bool show_shell,
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

                const double physical_core_radius = self.core_diameter / 2.0;
                const double physical_outer_radius = physical_core_radius + self.shell_thickness;

                double display_core_radius = 0.0;

                if (physical_outer_radius > 0.0) {
                    display_core_radius = display_radius * physical_core_radius / physical_outer_radius;
                }

                if (show_shell) {
                    add_sphere_surface_to_ax(
                        ax,
                        display_radius,
                        shell_color,
                        shell_alpha,
                        theta_sampling,
                        phi_sampling,
                        linewidth,
                        shade
                    );
                }

                if (show_core) {
                    add_sphere_surface_to_ax(
                        ax,
                        display_core_radius,
                        core_color,
                        core_alpha,
                        theta_sampling,
                        phi_sampling,
                        linewidth,
                        shade
                    );
                }

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

                set_coreshell_axis_limits(
                    ax,
                    axis_limit,
                    show_axes
                );
            },
            py::arg("ax"),
            py::arg("core_color") = py::str("black"),
            py::arg("shell_color") = py::str("grey"),
            py::arg("core_alpha") = 1.0,
            py::arg("shell_alpha") = 0.35,
            py::arg("display_radius") = 0.1,
            py::arg("show_core") = true,
            py::arg("show_shell") = true,
            py::arg("show_unit_sphere") = true,
            py::arg("unit_sphere_color") = py::str("grey"),
            py::arg("unit_sphere_alpha") = 0.15,
            py::arg("theta_sampling") = 24,
            py::arg("phi_sampling") = 24,
            py::arg("linewidth") = 0.0,
            py::arg("shade") = false,
            py::arg("show_axes") = false,
            R"pbdoc(
                Add the core shell geometry to a Matplotlib 3D axis.

                This private plotting helper replaces the previous PyVista based
                ``add_to_scene`` method. It draws the outer shell surface and,
                optionally, the core surface and a transparent unit sphere for
                angular reference.

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Matplotlib axis created with ``projection="3d"``.
                core_color : str or tuple, optional
                    Matplotlib color used for the core surface.
                shell_color : str or tuple, optional
                    Matplotlib color used for the shell surface.
                core_alpha : float, optional
                    Transparency of the core surface.
                shell_alpha : float, optional
                    Transparency of the shell surface.
                display_radius : float, optional
                    Outer radius used for visualization. This is a display
                    scale, not the physical radius of the particle.
                show_core : bool, optional
                    If ``True``, draw the core surface.
                show_shell : bool, optional
                    If ``True``, draw the shell outer surface.
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
                The core radius is scaled relative to ``display_radius`` using
                the physical ratio between the core radius and the outer particle
                radius. The absolute plotted size is therefore normalized for
                visualization, while the core shell proportions are preserved.
            )pbdoc"
        )
        ;
}