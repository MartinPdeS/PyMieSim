#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "./properties_set.h"

namespace py = pybind11;

void register_properties_set(py::module& module) {

    py::enum_<PropertyMode>(module, "PropertyMode")
        .value("Constant", PropertyMode::Constant)
        .value("Spectral", PropertyMode::Spectral)
        .export_values();

    py::class_<ScattererProperties, std::shared_ptr<ScattererProperties>>(module, "ScattererProperties")
        .def(py::init<std::vector<complex128>>(), py::arg("values"))
        .def(py::init<std::vector<std::vector<complex128>>>(), py::arg("values"))
        .def_property_readonly(
            "size",
            [](const ScattererProperties &self) {
                return self.size();
            }
        )
        .def("get", &ScattererProperties::get)
        .def("is_constant", &ScattererProperties::is_constant)
        .def("is_spectral", &ScattererProperties::is_spectral)
        .def(
            "__iter__",
            [](ScattererProperties &self) -> py::object {

                if (self.mode == PropertyMode::Constant) {
                    return py::iter(py::cast(self.constant_values));
                }

                return py::iter(py::cast(self.spectral_values));
            }
        )
        .def("__len__", [](const ScattererProperties &self) {
            return self.size();
        })
        .def(
            "__eq__",
            [](const ScattererProperties &a, const ScattererProperties &b) {

                if (a.mode != b.mode)
                    return false;

                if (a.mode == PropertyMode::Constant)
                    return a.constant_values == b.constant_values;

                return a.spectral_values == b.spectral_values;
            },
            py::is_operator()
        )
        ;


    py::class_<MediumProperties, std::shared_ptr<MediumProperties>>(module, "MediumProperties")
        .def(py::init<std::vector<double>>(), py::arg("values"))
        .def(py::init<std::vector<std::vector<double>>>(), py::arg("values"))
        .def_property_readonly(
            "size",
            [](const MediumProperties &self) {
                return self.size();
            }
        )
        .def("get", &MediumProperties::get)
        .def("is_constant", &MediumProperties::is_constant)
        .def("is_spectral", &MediumProperties::is_spectral)
        .def(
            "__iter__",
            [](MediumProperties &self) -> py::object {

                if (self.mode == PropertyMode::Constant) {
                    return py::iter(py::cast(self.constant_values));
                }

                return py::iter(py::cast(self.spectral_values));
            }
        )
        .def("__len__", [](const MediumProperties &self) {
            return self.size();
        })
        .def(
            "__eq__",
            [](const MediumProperties &a, const MediumProperties &b) {

                if (a.mode != b.mode)
                    return false;

                if (a.mode == PropertyMode::Constant)
                    return a.constant_values == b.constant_values;

                return a.spectral_values == b.spectral_values;
            },
            py::is_operator()
        )
        ;
}