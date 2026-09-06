#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <sstream>
#include <pybind11/complex.h>

#include <pint/pint.h>
#include "./material_set.h"

namespace py = pybind11;


PYBIND11_MODULE(material_set, module) {
    py::object ureg = get_shared_ureg();

    auto format_material_like_repr = [](const std::string& name, size_t size) {
        std::ostringstream stream;
        stream << "<" << name << " elements=" << size << ">";
        return stream.str();
    };

    py::class_<MaterialSet, std::shared_ptr<MaterialSet>>(
        module,
        "MaterialSet",
        R"pdoc(
            Container of material definitions used by experiment scatterer and
            detector sets.

            A ``MaterialSet`` behaves like a lightweight sequence of material
            objects or refractive index values converted to material objects.
        )pdoc"
    )
        .def(py::init<>())
        .def(
            py::init<const std::vector<std::shared_ptr<BaseMaterial>>&>(),
            py::arg("materials"),
            R"pdoc(
                Initialize a material set from existing material objects.

                Parameters
                ----------
                materials : list[BaseMaterial]
                    Material instances stored in the set.
            )pdoc"
        )
        .def(
            py::init<const std::vector<std::complex<double>>&>(),
            py::arg("refractive_indices"),
            R"pdoc(
                Initialize a material set from complex refractive index values.

                Parameters
                ----------
                refractive_indices : list[complex]
                    Complex refractive indices converted into constant material
                    objects.
            )pdoc"
        )
        .def("__len__", [](const MaterialSet& self){return self.size();}, "Return the number of stored materials.")
        .def("__getitem__", [](const MaterialSet& self, std::size_t index) {
            if (index >= self.size()) {
                throw py::index_error();
            }
            return self[index];
        }, py::arg("index"), "Return the material at the given index.")
        .def(
            "__repr__",
            [format_material_like_repr](const MaterialSet& self) {
                return format_material_like_repr("MaterialSet", self.size());
            },
            "Return a compact representation showing the number of stored materials."
        )
        ;

    py::class_<MediumSet, std::shared_ptr<MediumSet>>(
        module,
        "MediumSet",
        R"pdoc(
            Container of medium definitions used by experiment scatterer and
            detector sets.

            A ``MediumSet`` behaves like a lightweight sequence of medium
            objects or refractive index values converted to medium objects.
        )pdoc"
    )
        .def(py::init<>())
        .def(
            py::init<const std::vector<std::shared_ptr<BaseMedium>>&>(),
            py::arg("mediums"),
            R"pdoc(
                Initialize a medium set from existing medium objects.

                Parameters
                ----------
                mediums : list[BaseMedium]
                    Medium instances stored in the set.
            )pdoc"
        )
        .def(
            py::init<const std::vector<double>&>(),
            py::arg("refractive_indices"),
            R"pdoc(
                Initialize a medium set from refractive index values.

                Parameters
                ----------
                refractive_indices : list[float]
                    Refractive indices converted into constant medium objects.
            )pdoc"
        )
        .def("__len__", [](const MediumSet& self){return self.size();}, "Return the number of stored media.")
        .def("__getitem__", [](const MediumSet& self, std::size_t index) {
            if (index >= self.size()) {
                throw py::index_error();
            }
            return self[index];
        }, py::arg("index"), "Return the medium at the given index.")
        .def(
            "__repr__",
            [format_material_like_repr](const MediumSet& self) {
                return format_material_like_repr("MediumSet", self.size());
            },
            "Return a compact representation showing the number of stored media."
        )
        ;
}

// -
