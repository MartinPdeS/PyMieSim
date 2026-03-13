#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "./base.h"
#include "./constant.h"
#include "./tabulated.h"
#include "./sellmeier.h"
#include <pint/pint.h>

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace py = pybind11;


PYBIND11_MODULE(material, module) {
    py::object ureg = get_shared_ureg();

    module.def(
        "print_available",
        []() {
            py::module_ pyoptik = py::module_::import("PyOptik");
            py::object Material = pyoptik.attr("Material");

            Material.attr("use_tabulated") = py::bool_(true);
            Material.attr("use_sellmeier") = py::bool_(true);

            Material.attr("print_available")();
        }
    );

    py::class_<BaseMaterial, std::shared_ptr<BaseMaterial>>(module, "BaseMaterial")
        .def(
            "initialize",
            [](BaseMaterial& self, const py::object& wavelength) {
                return self.initialize(
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>()
                );
            }
        )
        .def(
            "get_refractive_index",
            [ureg](const BaseMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        ;

    py::class_<ConstantMaterial, BaseMaterial, std::shared_ptr<ConstantMaterial>>(module, "ConstantMaterial")
        .def(
            py::init([ureg](const py::object& refractive_index) {
                const complex128 refractive_index_value =
                    refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<complex128>();

                return std::make_shared<ConstantMaterial>(refractive_index_value);
            }),
            py::arg("refractive_index")
        )
        .def_property_readonly(
            "refractive_index",
            [ureg](const ConstantMaterial& self) {
                py::object magnitude = py::cast(self.constant_refractive_index);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def(
            "get_refractive_index",
            [ureg](const ConstantMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def(
            "__repr__",
            [](const ConstantMaterial& self) {
                char buffer[128];

                if (self.constant_refractive_index.imag() >= 0.0) {
                    std::snprintf(
                        buffer,
                        sizeof(buffer),
                        "%g + i%g",
                        self.constant_refractive_index.real(),
                        self.constant_refractive_index.imag()
                    );
                }
                else {
                    std::snprintf(
                        buffer,
                        sizeof(buffer),
                        "%g - i%g",
                        self.constant_refractive_index.real(),
                        -self.constant_refractive_index.imag()
                    );
                }

                return std::string(buffer);
            }
        )
        ;

    py::class_<TabulatedMaterial, BaseMaterial, std::shared_ptr<TabulatedMaterial>>(module, "TabulatedMaterial")
        .def(
            py::init([](const py::str& material_name) {
                py::module_ pyoptik = py::module_::import("PyOptik");
                py::object Material = pyoptik.attr("Material");

                Material.attr("use_tabulated") = py::bool_(true);
                Material.attr("use_sellmeier") = py::bool_(false);

                py::object material = Material.attr(material_name);

                std::vector<double> n_values = material.attr("n_values").cast<std::vector<double>>();
                std::vector<double> k_values = material.attr("k_values").cast<std::vector<double>>();

                std::vector<complex128> refractive_index_values;
                refractive_index_values.reserve(n_values.size());

                for (std::size_t index = 0; index < n_values.size(); ++index) {
                    refractive_index_values.emplace_back(n_values[index], k_values[index]);
                }

                return std::make_shared<TabulatedMaterial>(
                    material_name.cast<std::string>(),
                    material.attr("wavelength").attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                    refractive_index_values,
                    true
                );
            }),
            py::arg("material_name")
        )
        .def(
            py::init([ureg](
                const std::string& name,
                const py::object& wavelengths,
                const py::object& refractive_indices,
                const bool allow_extrapolation
            ) {
                const std::vector<double> wavelength_values =
                    wavelengths.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                const std::vector<complex128> refractive_index_values =
                    refractive_indices.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<std::vector<complex128>>();

                return std::make_shared<TabulatedMaterial>(
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
        .def_readonly("name", &TabulatedMaterial::name)
        .def_property_readonly(
            "wavelengths",
            [ureg](const TabulatedMaterial& self) {
                py::object magnitude = py::cast(self.wavelengths);
                return magnitude * ureg.attr("meter");
            }
        )
        .def_property_readonly(
            "refractive_indices",
            [ureg](const TabulatedMaterial& self) {
                py::object magnitude = py::cast(self.refractive_indices);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def_readonly("allow_extrapolation", &TabulatedMaterial::allow_extrapolation)
        .def(
            "get_refractive_index",
            [ureg](const TabulatedMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def(
            "__repr__",
            [](const TabulatedMaterial& self) {
                return self.name;
            }
        )
        ;

    py::class_<SellmeierMaterial, BaseMaterial, std::shared_ptr<SellmeierMaterial>>(module, "SellmeierMaterial")
        .def(
            py::init([](const py::str& material_name) {
                py::module_ pyoptik = py::module_::import("PyOptik");
                py::object Material = pyoptik.attr("Material");

                Material.attr("use_tabulated") = py::bool_(false);
                Material.attr("use_sellmeier") = py::bool_(true);

                py::object material = Material.attr(material_name);

                return std::make_shared<SellmeierMaterial>(
                    material_name.cast<std::string>(),
                    material.attr("coefficients").cast<std::vector<double>>(),
                    material.attr("formula_type").cast<std::size_t>(),
                    material.attr("wavelength_bound").attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                    true
                );
            }),
            py::arg("material_name")
        )
        .def(
            py::init([ureg](
                const std::string& name,
                const std::vector<double>& coefficients,
                const std::size_t formula_type,
                const py::object& wavelength_bound,
                const bool allow_extrapolation
            ) {
                const std::vector<double> wavelength_bound_values =
                    wavelength_bound.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                return std::make_shared<SellmeierMaterial>(
                    name,
                    coefficients,
                    formula_type,
                    wavelength_bound_values,
                    allow_extrapolation
                );
            }),
            py::arg("name"),
            py::arg("coefficients"),
            py::arg("formula_type"),
            py::arg("wavelength_bound") = py::cast(std::vector<double>{}) * ureg.attr("meter"),
            py::arg("allow_extrapolation") = false
        )
        .def_readonly("name", &SellmeierMaterial::name)
        .def_readonly("coefficients", &SellmeierMaterial::coefficients)
        .def_readonly("formula_type", &SellmeierMaterial::formula_type)
        .def_readonly("allow_extrapolation", &SellmeierMaterial::allow_extrapolation)
        .def_property_readonly(
            "wavelength_bound",
            [ureg](const SellmeierMaterial& self) {
                py::object magnitude = py::cast(self.wavelength_bound);
                return magnitude * ureg.attr("meter");
            }
        )
        .def(
            "get_refractive_index",
            [ureg](const SellmeierMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def(
            "__repr__",
            [](const SellmeierMaterial& self) {
                return self.name;
            }
        )
        ;

    py::class_<BaseMedium, std::shared_ptr<BaseMedium>>(module, "BaseMedium")
        .def(
            "get_refractive_index",
            [ureg](const BaseMedium& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        ;

    py::class_<ConstantMedium, BaseMedium, std::shared_ptr<ConstantMedium>>(module, "ConstantMedium")
        .def(
            py::init([ureg](const py::object& refractive_index) {
                const double refractive_index_value =
                    refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();

                return std::make_shared<ConstantMedium>(refractive_index_value);
            }),
            py::arg("refractive_index")
        )
        .def_property_readonly(
            "refractive_index",
            [ureg](const ConstantMedium& self) {
                py::object magnitude = py::cast(self.constant_refractive_index);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def(
            "get_refractive_index",
            [ureg](const ConstantMedium& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def(
            "__repr__",
            [](const ConstantMedium& self) {
                char buffer[64];
                std::snprintf(buffer, sizeof(buffer), "%g", self.constant_refractive_index);
                return std::string(buffer);
            }
        )
        ;

    py::class_<TabulatedMedium, BaseMedium, std::shared_ptr<TabulatedMedium>>(module, "TabulatedMedium")
        .def(
            py::init([ureg](
                const std::string& name,
                const py::object& wavelengths,
                const py::object& refractive_indices,
                const bool allow_extrapolation
            ) {
                const std::vector<double> wavelength_values =
                    wavelengths.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                const std::vector<double> refractive_index_values =
                    refractive_indices.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<std::vector<double>>();

                return std::make_shared<TabulatedMedium>(
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
        .def_readonly("name", &TabulatedMedium::name)
        .def_property_readonly(
            "wavelengths",
            [ureg](const TabulatedMedium& self) {
                py::object magnitude = py::cast(self.wavelengths);
                return magnitude * ureg.attr("meter");
            }
        )
        .def_property_readonly(
            "refractive_indices",
            [ureg](const TabulatedMedium& self) {
                py::object magnitude = py::cast(self.refractive_indices);
                return magnitude * ureg.attr("RIU");
            }
        )
        .def_readonly("allow_extrapolation", &TabulatedMedium::allow_extrapolation)
        .def(
            "get_refractive_index",
            [ureg](const TabulatedMedium& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def(
            "__repr__",
            [](const TabulatedMedium& self) {
                return self.name;
            }
        )
        ;

    py::class_<SellmeierMedium, BaseMedium, std::shared_ptr<SellmeierMedium>>(module, "SellmeierMedium")
        .def(
            py::init([](const py::str& material_name) {
                py::module_ pyoptik = py::module_::import("PyOptik");
                py::object Material = pyoptik.attr("Material");

                Material.attr("use_tabulated") = py::bool_(false);
                Material.attr("use_sellmeier") = py::bool_(true);

                py::object material = Material.attr(material_name);

                return std::make_shared<SellmeierMedium>(
                    material_name.cast<std::string>(),
                    material.attr("coefficients").cast<std::vector<double>>(),
                    material.attr("formula_type").cast<std::size_t>(),
                    material.attr("wavelength_bound").attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                    true
                );
            }),
            py::arg("material_name")
        )
        .def(
            py::init([ureg](
                const std::string& name,
                const std::vector<double>& coefficients,
                const std::size_t formula_type,
                const py::object& wavelength_bound,
                const bool allow_extrapolation
            ) {
                const std::vector<double> wavelength_bound_values =
                    wavelength_bound.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                return std::make_shared<SellmeierMedium>(
                    name,
                    coefficients,
                    formula_type,
                    wavelength_bound_values,
                    allow_extrapolation
                );
            }),
            py::arg("name"),
            py::arg("coefficients"),
            py::arg("formula_type"),
            py::arg("wavelength_bound") = py::cast(std::vector<double>{}) * ureg.attr("meter"),
            py::arg("allow_extrapolation") = false
        )
        .def_readonly("name", &SellmeierMedium::name)
        .def_readonly("coefficients", &SellmeierMedium::coefficients)
        .def_readonly("formula_type", &SellmeierMedium::formula_type)
        .def_readonly("allow_extrapolation", &SellmeierMedium::allow_extrapolation)
        .def_property_readonly(
            "wavelength_bound",
            [ureg](const SellmeierMedium& self) {
                py::object magnitude = py::cast(self.wavelength_bound);
                return magnitude * ureg.attr("meter");
            }
        )
        .def(
            "get_refractive_index",
            [ureg](const SellmeierMedium& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude * ureg.attr("RIU");
            },
            py::arg("wavelength")
        )
        .def(
            "__repr__",
            [](const SellmeierMedium& self) {
                return self.name;
            }
        )
        ;
}