#include <pybind11/pybind11.h>
#include "./coreshell.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>

namespace py = pybind11;

void register_coreshell(py::module_& module) {
    py::object ureg = get_shared_ureg();

    // Binding for CoreShell class
    py::class_<CoreShell, BaseScatterer, std::shared_ptr<CoreShell>>(module, "CoreShell")
        .def(
            py::init([ureg](
                py::object core_diameter,
                py::object shell_thickness,
                py::object core_refractive_index,
                py::object shell_refractive_index,
                py::object medium_refractive_index,
                std::shared_ptr<BaseSource> source
            ) {
                py::object units_length = py::module_::import("PyMieSim.units").attr("Length");
                py::object units_riu = py::module_::import("PyMieSim.units").attr("RefractiveIndex");

                core_diameter = units_length.attr("check")(core_diameter);
                shell_thickness = units_length.attr("check")(shell_thickness);
                core_refractive_index = units_riu.attr("check")(core_refractive_index);
                shell_refractive_index = units_riu.attr("check")(shell_refractive_index);
                medium_refractive_index = units_riu.attr("check")(medium_refractive_index);
                double core_diameter_meter =
                    core_diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                double shell_thickness_meter =
                    shell_thickness.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                std::complex<double> core_refractive_index_riu =
                    core_refractive_index.attr("to")(ureg.attr("RIU")) .attr("magnitude") .cast<std::complex<double>>();

                std::complex<double> shell_refractive_index_riu =
                    shell_refractive_index.attr("to")(ureg.attr("RIU")) .attr("magnitude") .cast<std::complex<double>>();

                double medium_refractive_index_value =
                    medium_refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();

                return std::make_shared<CoreShell>(
                    core_diameter_meter,
                    shell_thickness_meter,
                    core_refractive_index_riu,
                    shell_refractive_index_riu,
                    medium_refractive_index_value,
                    std::move(source)
                );
            }),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            R"pbdoc(
                Constructor for CORESHELL, initializing it with physical and optical properties.

                Parameters
                ----------
                core_diameter : float
                    The diameter of the core of the spherical shell.
                shell_thickness : float
                    The thickness of the shell surrounding the core.
                core_refractive_index : complex
                    The refractive index of the core.
                shell_refractive_index : complex
                    The refractive index of the shell.
                medium_refractive_index : float
                    The refractive index of the surrounding medium.
                source : BaseSource
                    The source of the incident light.
            )pbdoc"
        )
        .def_readonly(
            "source",
            &CoreShell::source,
            "Source of the core-shell scatterer."
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
        .def_property_readonly(
            "core_refractive_index",
            [ureg](const CoreShell &self) {
                py::object magnitude = py::cast(self.core_refractive_index);
                return (magnitude * ureg.attr("RIU")).attr("to_compact")();
            },
            "Refractive index of the core shell."
        )
        .def_property_readonly(
            "shell_refractive_index",
            [ureg](const CoreShell &self) {
                py::object magnitude = py::cast(self.shell_refractive_index);
                return (magnitude * ureg.attr("RIU")).attr("to_compact")();
            },
            "Refractive index of the shell."
        )
        .def_readonly(
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
        ;
}
