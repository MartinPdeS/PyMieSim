#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "./material.h"
#include "./medium.h"
#include <pint/pint.h>

namespace py = pybind11;

PYBIND11_MODULE(material, module) {
    py::object ureg = get_shared_ureg();

    py::class_<BaseMaterial, std::shared_ptr<BaseMaterial>>(module, "BaseMaterial")
        .def(
            "get_refractive_index",
            [ureg](const BaseMaterial& self, py::object wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<double>();

                py::object magnitude = py::cast(
                    self.get_refractive_index(wavelength_meter)
                );

                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        );

    py::class_<ConstantMaterial, BaseMaterial, std::shared_ptr<ConstantMaterial>>(module, "ConstantMaterial")
        .def(
            py::init([ureg](py::object refractive_index) {
                const complex128 refractive_index_value =
                    refractive_index.attr("to")(ureg.attr("RIU"))
                    .attr("magnitude")
                    .cast<complex128>();

                return std::make_shared<ConstantMaterial>(
                    refractive_index_value
                );
            }),
            py::arg("refractive_index")
        )
        .def_property_readonly(
            "refractive_index",
            [ureg](const ConstantMaterial& self) {
                py::object magnitude = py::cast(self.refractive_index);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def(
            "get_refractive_index",
            [ureg](const ConstantMaterial& self, py::object wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<double>();

                py::object magnitude = py::cast(
                    self.get_refractive_index(wavelength_meter)
                );

                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def("__repr__", [](const ConstantMaterial& self) {
            return std::to_string(self.constant_refractive_index.real()) + " + " + std::to_string(self.constant_refractive_index.imag()) + "j>";
        })
        ;

    py::class_<Material, BaseMaterial, std::shared_ptr<Material>>(module, "Material")
        .def(
            py::init([ureg](
                const std::string& name,
                py::object wavelengths,
                py::object refractive_indices,
                bool allow_extrapolation
            ) {
                const std::vector<double> wavelength_values =
                    wavelengths.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                const std::vector<complex128> refractive_index_values =
                    refractive_indices.attr("to")(ureg.attr("RIU"))
                    .attr("magnitude")
                    .cast<std::vector<complex128>>();

                return std::make_shared<Material>(
                    name,
                    wavelength_values,
                    refractive_index_values,
                    allow_extrapolation
                );
            }),
            py::arg("name"),
            py::arg("wavelengths"),
            py::arg("refractive_indices"),
            py::arg("allow_extrapolation") = false
        )
        .def_readonly(
            "name",
            &Material::name
        )
        .def_property_readonly(
            "refractive_indices",
            [ureg](const Material& self) {
                py::object magnitude = py::cast(self.refractive_indices);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def_property_readonly(
            "wavelengths",
            [ureg](const Material& self) {
                py::object magnitude = py::cast(self.wavelengths);
                return magnitude * ureg.attr("meter");
            }
        )
        .def_readonly(
            "allow_extrapolation",
            &Material::allow_extrapolation
        )
        .def(
            "get_refractive_index",
            [ureg](const Material& self, py::object wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<double>();

                py::object magnitude = py::cast(
                    self.get_refractive_index(wavelength_meter)
                );

                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def("__repr__", [](const Material& self) {
            return self.name;
        }
    )
    ;

    py::class_<BaseMedium, std::shared_ptr<BaseMedium>>(module, "BaseMedium")
        .def(
            "get_refractive_index",
            [ureg](const BaseMedium& self, py::object wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<double>();

                py::object magnitude = py::cast(
                    self.get_refractive_index(wavelength_meter)
                );

                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        );

    py::class_<ConstantMedium, BaseMedium, std::shared_ptr<ConstantMedium>>(module, "ConstantMedium")
        .def(
            py::init([ureg](py::object refractive_index) {
                const double refractive_index_value =
                    refractive_index.attr("to")(ureg.attr("RIU"))
                    .attr("magnitude")
                    .cast<double>();

                return std::make_shared<ConstantMedium>(
                    refractive_index_value
                );
            }),
            py::arg("refractive_index")
        )
        .def_property_readonly(
            "refractive_index",
            [ureg](const ConstantMedium& self) {
                py::object magnitude = py::cast(self.refractive_index);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def(
            "get_refractive_index",
            [ureg](const ConstantMedium& self, py::object wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<double>();

                py::object magnitude = py::cast(
                    self.get_refractive_index(wavelength_meter)
                );

                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def("__repr__", [](const ConstantMedium& self) {
            return std::to_string(self.constant_refractive_index);
        })
        ;

    py::class_<Medium, BaseMedium, std::shared_ptr<Medium>>(module, "Medium")
        .def(
            py::init([ureg](
                const std::string& name,
                py::object wavelengths,
                py::object refractive_indices,
                bool allow_extrapolation
            ) {
                const std::vector<double> wavelength_values =
                    wavelengths.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                const std::vector<double> refractive_index_values =
                    refractive_indices.attr("to")(ureg.attr("RIU"))
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                return std::make_shared<Medium>(
                    name,
                    wavelength_values,
                    refractive_index_values,
                    allow_extrapolation
                );
            }),
            py::arg("name"),
            py::arg("wavelengths"),
            py::arg("refractive_indices"),
            py::arg("allow_extrapolation") = false
        )
        .def_readonly(
            "name",
            &Medium::name
        )
        .def_property_readonly(
            "refractive_indices",
            [ureg](const Medium& self) {
                py::object magnitude = py::cast(self.refractive_indices);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def_property_readonly(
            "wavelengths",
            [ureg](const Medium& self) {
                py::object magnitude = py::cast(self.wavelengths);
                return magnitude * ureg.attr("meter");
            }
        )
        .def_readonly(
            "allow_extrapolation",
            &Medium::allow_extrapolation
        )
        .def(
            "get_refractive_index",
            [ureg](const Medium& self, py::object wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter"))
                    .attr("magnitude")
                    .cast<double>();

                py::object magnitude = py::cast(
                    self.get_refractive_index(wavelength_meter)
                );

                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def("__repr__", [](const Medium& self) {
            return self.name;
        })
        ;
}