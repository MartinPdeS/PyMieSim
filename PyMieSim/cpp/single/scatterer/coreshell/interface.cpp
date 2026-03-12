#include <pybind11/pybind11.h>
#include "./coreshell.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>
#include <single/scatterer/utils.h>

namespace py = pybind11;

void register_coreshell(py::module_& module) {
    py::object ureg = get_shared_ureg();

    py::class_<CoreShell, BaseScatterer, std::shared_ptr<CoreShell>>(module, "CoreShell")
        .def(
            py::init([ureg](
                py::object core_diameter,
                py::object shell_thickness,
                const py::object& core_material,
                const py::object& shell_material,
                const py::object& medium,
                int max_order = 0
            ) {
                double core_diameter_meter =
                    core_diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                double shell_thickness_meter =
                    shell_thickness.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

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
            }),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_material"),
            py::arg("shell_material"),
            py::arg("medium"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Constructor for CORESHELL, initializing it with physical and optical properties.

                Parameters
                ----------
                core_diameter : float
                    The diameter of the core of the spherical shell.
                shell_thickness : float
                    The thickness of the shell surrounding the core.
                core_material : BaseMaterial
                    The material of the core.
                shell_material : BaseMaterial
                    The material of the shell.
                medium : BaseMedium
                    The surrounding medium of the core-shell scatterer.
                max_order : int, optional
                    The maximum order of the scattering coefficients to compute. If set to 0, it will
            )pbdoc"
        )
        .def(
            "print_properties",
            &CoreShell::print_properties,
            "Prints the properties of the core-shell scatterer."
        )
        .def_property_readonly(
            "core_diameter",
            [ureg](const CoreShell &self) {
                return (py::float_(self.core_diameter) * ureg.attr("meter")).attr("to_compact")();
            },
            "Diameter of the core shell."
        )
        .def_property_readonly(
            "shell_thickness",
            [ureg](const CoreShell &self) {
                return (py::float_(self.shell_thickness) * ureg.attr("meter")).attr("to_compact")();
            },
            "Thickness of the shell."
        )
        .def_readonly(
            "core_material",
            &CoreShell::core_material,
            "Material of the core shell."
        )
        .def_readonly(
            "shell_material",
            &CoreShell::shell_material,
            "Material of the shell."
        )
        .def_readonly_static(
            "property_names",
            &CoreShell::property_names,
            "Property names of the core-shell scatterer."
        )
        .def_property_readonly(
            "radius",
            [ureg](const CoreShell &self) {
                double radius = self.core_diameter / 2.0 + self.shell_thickness;
                return (py::float_(radius) * ureg.attr("meter")).attr("to_compact")();
            },
            "Overall radius of the core-shell scatterer."
        )
        .def_property_readonly(
            "volume",
            [ureg](const CoreShell &self) {
                double core_radius = self.core_diameter / 2.0;
                double shell_outer_radius = core_radius + self.shell_thickness;
                double volume = (4.0 / 3.0) * Constants::PI * (shell_outer_radius * shell_outer_radius * shell_outer_radius - core_radius * core_radius * core_radius);
                return (py::float_(volume) * ureg.attr("meter**3")).attr("to_compact")();
            },
            "Volume of the core-shell scatterer."
        )
        .def_property("an",
            [ureg](CoreShell& self) {return vector_as_numpy_view(self, self.an);},
            [ureg](CoreShell& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.an, arr);
            },
            R"pbdoc(
                Returns the 'an' scattering coefficients.

                .. math::
                    a_n = \frac{
                        \Psi_n (y)
                        \left[ \dot{\Psi}_n (\tilde{m}_{2} y) - A_n (x) \dot{\xi}_n (\tilde{m}_{2} y) \right]
                        -
                        \tilde{m}_{2} \dot{\Psi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - A_n (x) \xi_n (\tilde{m}_{2} y) \right]
                    }{
                        \xi_n (y)
                        \left[ \dot{\Psi}_n (\tilde{m}_{2} y) - A_n (x) \dot{\xi}_n (\tilde{m}_{2} y) \right]
                        -
                        \tilde{m}_{2} \dot{\xi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - A_n (x) \xi_n (\tilde{m}_{2} y) \right]
                    }

                with:

                .. math::
                    A_n = \frac{
                        \tilde{m}_{2}  \Psi_n       (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\Psi}_n (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \xi_n        (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\xi}_n  (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }
                    \\
                    B_n = \frac{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\psi}_n (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \psi_n  (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\xi}_n  (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \xi_n   (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }

                Returns
                -------
                list
                    A list of 'an' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def_property("bn",
            [ureg](CoreShell& self) {return vector_as_numpy_view(self, self.bn);},
            [ureg](CoreShell& self,
                py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                    vector_assign_from_numpy(self.bn, arr);
            },
            R"pbdoc(
                Returns the 'bn' scattering coefficients.

                .. math::
                    b_n = \frac{
                        \tilde{m}_{2}
                        \Psi_n (y)
                        \left[ \dot{\Psi}_n  (\tilde{m}_{2} y) - B_n (x) \dot{\xi}_n  (\tilde{m}_{2} y) \right]
                        -
                        \dot{\Psi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - B_n (x) \xi_n  (\tilde{m}_{2} y) \right]
                    }{
                        \tilde{m}_{2}
                        \xi_n (y)
                        \left[ \dot{\Psi}_n  (\tilde{m}_{2} y) - A_n (x) \dot{\xi}_n  (\tilde{m}_{2} y) \right] -
                        \dot{\xi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - B_n (x) \xi_n (\tilde{m}_{2} y) \right]
                    }

                with:

                .. math::
                    A_n = \frac{
                        \tilde{m}_{2}  \Psi_n       (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\Psi}_n (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \xi_n        (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\xi}_n  (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }
                    \\
                    B_n = \frac{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\psi}_n (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \psi_n  (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\xi}_n  (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \xi_n   (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }

                Returns
                -------
                list
                    A list of 'bn' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def_property("cn",
            [ureg](CoreShell& self) {return vector_as_numpy_view(self, self.cn);},
            [ureg](CoreShell& self,
                py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                    vector_assign_from_numpy(self.cn, arr);
            },
            R"pbdoc(
                Returns the 'cn' scattering coefficients.

                Returns
                -------
                list
                    A list of 'cn' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def_property("dn",
            [ureg](CoreShell& self) {return vector_as_numpy_view(self, self.dn);},
            [ureg](CoreShell& self,
                py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                    vector_assign_from_numpy(self.dn, arr);
            },
            R"pbdoc(
                Returns the 'dn' scattering coefficients.

                Returns
                -------
                list
                    A list of 'dn' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def(
            "add_to_scene",
            [](const CoreShell& self, py::object& scene, const py::kwargs& kwargs = py::dict()) {
                py::module_ pv = py::module_::import("pyvista");
                py::object sphere = pv.attr("Sphere");


                py::kwargs sphere_kwargs = py::kwargs();
                sphere_kwargs["radius"] = self.core_diameter / 2.0 + self.shell_thickness;
                sphere_kwargs["center"] = py::make_tuple(0, 0, 0);
                sphere_kwargs["theta_resolution"] = 50;
                sphere_kwargs["phi_resolution"] = 50;

                py::object py_sphere = sphere(**sphere_kwargs);
                scene.attr("add_mesh")(py_sphere, **kwargs);

            }
        )
        .def(
            "add_to_scene",
            [](const CoreShell& self, const py::object& scene, py::kwargs kwargs) {
                py::module_ pyvista = py::module_::import("pyvista");

                py::object physical_sphere = pyvista.attr("Sphere")(
                    py::arg("radius") = py::float_(0.1),
                    py::arg("center") = py::make_tuple(0.0, 0.0, 0.0),
                    py::arg("theta_resolution") = 50,
                    py::arg("phi_resolution") = 50
                );

                py::object unit_sphere = pyvista.attr("Sphere")(
                    py::arg("radius") = py::float_(1.0),
                    py::arg("center") = py::make_tuple(0.0, 0.0, 0.0),
                    py::arg("theta_resolution") = 50,
                    py::arg("phi_resolution") = 50
                );

                scene.attr("add_mesh")(physical_sphere, **kwargs);

                py::dict unit_sphere_kwargs;
                unit_sphere_kwargs["opacity"] = py::float_(0.15);

                if (kwargs.contains("color")) {
                    unit_sphere_kwargs["color"] = kwargs["color"];
                }
                else {
                    unit_sphere_kwargs["color"] = py::str("white");
                }

                scene.attr("add_mesh")(unit_sphere, **unit_sphere_kwargs);
            },
            py::arg("scene")
        )
        ;
}
