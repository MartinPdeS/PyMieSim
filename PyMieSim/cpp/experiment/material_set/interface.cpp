#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers

#include <pint/pint.h>
#include "./material_set.h"
#include "./medium_set.h"

namespace py = pybind11;


PYBIND11_MODULE(material_set, module) {
    py::object ureg = get_shared_ureg();

    py::class_<MaterialSet, std::shared_ptr<MaterialSet>>(module, "MaterialSet")
        .def(py::init<>())
        .def(py::init<const std::vector<std::shared_ptr<BaseMaterial>>&>())
        .def(py::init<const std::vector<std::complex<double>>&>())
        .def("__len__", [](const MaterialSet& self){return self.size();})
        .def("__getitem__", [](const MaterialSet& self, std::size_t index) {
            if (index >= self.size()) {
                throw py::index_error();
            }
            return self[index];
        })
        ;

    py::class_<MediumSet, std::shared_ptr<MediumSet>>(module, "MediumSet")
        .def(py::init<>())
        .def(py::init<const std::vector<std::shared_ptr<BaseMedium>>&>())
        .def(py::init<const std::vector<double>&>())
        .def("__len__", [](const MediumSet& self){return self.size();})
        .def("__getitem__", [](const MediumSet& self, std::size_t index) {
            if (index >= self.size()) {
                throw py::index_error();
            }
            return self[index];
        })
        ;
}

// -
