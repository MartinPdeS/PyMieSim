#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h>

namespace py = pybind11;


template <typename MaterialSetType, typename RefractiveIndexType, typename BaseClass>
MaterialSetType create_material_set_from_pyobject(
    const py::object& material_object,
    const std::string& material_name
) {
    if (py::isinstance<MaterialSetType>(material_object)) {
        return py::cast<MaterialSetType>(material_object);
    }

    if (py::isinstance<BaseClass>(material_object)) {
        return MaterialSetType(
            std::vector<std::shared_ptr<BaseClass>>{
                py::cast<std::shared_ptr<BaseClass>>(material_object)
            }
        );
    }

    if (
        py::isinstance<py::float_>(material_object) ||
        py::isinstance<py::int_>(material_object) ||
        py::isinstance<complex128>(material_object)
    ) {
        return MaterialSetType(
            std::vector<RefractiveIndexType>{
                material_object.cast<RefractiveIndexType>()
            }
        );
    }

    if (py::isinstance<py::sequence>(material_object) && !py::isinstance<py::str>(material_object)) {
        try {
            return MaterialSetType(
                py::cast<std::vector<RefractiveIndexType>>(material_object)
            );
        }
        catch (const py::cast_error&) {
        }

        try {
            return MaterialSetType(
                py::cast<std::vector<std::shared_ptr<BaseClass>>>(material_object)
            );
        }
        catch (const py::cast_error&) {
        }

        throw std::runtime_error(
            "Invalid type for " + material_name + ". Expected a MaterialSet, "
            "a scalar refractive index, a sequence of refractive indices, "
            "a Material instance, or a sequence of Material instances."
        );
    }

    throw std::runtime_error(
        "Invalid type for " + material_name + ". Expected a MaterialSet, "
        "a scalar refractive index, a sequence of refractive indices, "
        "a Material instance, or a sequence of Material instances."
    );
}