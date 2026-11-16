#include "pint.h"

#include <stdexcept>

namespace py = pybind11;

py::object get_shared_ureg() {
    py::gil_scoped_acquire gil;
    return py::module_::import("PyMieSim._pint").attr("get_ureg")();
}

py::object validated_quantity(py::handle value, const std::string &unit, const std::string &name, bool scalar) {
    py::object quantity = py::reinterpret_borrow<py::object>(value);
    py::module_::import("TypedUnit")
        .attr("validate_quantity")(quantity, py::arg("unit") = unit, py::arg("name") = name,
                                   py::arg("registry") = get_shared_ureg(), py::arg("scalar") = scalar);
    return quantity;
}

py::object validated_no_unit(py::handle value, const std::string &name, bool scalar) {
    py::module_ typedunit = py::module_::import("TypedUnit");
    py::object markers = scalar ? py::object(py::make_tuple(typedunit.attr("Scalar"), typedunit.attr("Finite")))
                                : py::object(py::make_tuple(typedunit.attr("Finite")));
    py::object annotation = typedunit.attr("NoUnit").attr("__class_getitem__")(markers);
    typedunit.attr("validate_annotation")(value, py::arg("annotation") = annotation, py::arg("name") = name);
    return py::reinterpret_borrow<py::object>(value);
}

double to_meters_strict(py::handle value) {
    py::object q = validated_quantity(value, "meter", "length", true);
    py::object q_m = q.attr("to")("meter");

    return py::cast<double>(q_m.attr("magnitude"));
}

py::object registry_from_quantity(py::handle value) {
    // Pint Quantity objects expose a registry as _REGISTRY
    if (!py::hasattr(value, "_REGISTRY")) {
        throw py::type_error("length must be a pint.Quantity (missing _REGISTRY)");
    }
    return py::reinterpret_borrow<py::object>(value).attr("_REGISTRY");
}

py::object meters_quantity_with_ureg(const py::object& ureg, double meters_value) {
    return ureg.attr("Quantity")(meters_value, "meter");
}


std::vector<double> to_vector_units(py::handle values, const std::string& unit) {
    if (!py::hasattr(values, "__iter__")) {
        throw py::type_error("expected an iterable of pint.Quantity objects");
    }

    std::vector<double> magnitudes;

    // Reserve only when __len__ exists to avoid surprising errors for generators
    if (py::hasattr(values, "__len__")) {
        magnitudes.reserve(py::len(values));
    }

    for (py::handle item : py::reinterpret_borrow<py::iterable>(values)) {
        py::object quantity = validated_quantity(item, unit, "values", true);
        py::object converted = quantity.attr("to")(unit);

        magnitudes.push_back(py::cast<double>(converted.attr("magnitude")));
    }

    return magnitudes;
}


double quantity_scalar_to_meters(py::object quantity) {
    py::object ureg = get_shared_ureg();
    py::object meter = ureg.attr("meter");
    validated_quantity(quantity, "meter", "length", true);
    py::object converted = quantity.attr("to")(meter);
    return py::float_(converted.attr("magnitude"));
}

std::vector<double> quantity_1d_to_meters_vector(py::object quantity) {
    py::object ureg = get_shared_ureg();
    py::object meter = ureg.attr("meter");

    validated_quantity(quantity, "meter", "radii");
    py::object converted = quantity.attr("to")(meter);
    py::object magnitude = converted.attr("magnitude");

    py::module_ numpy = py::module_::import("numpy");
    py::object float64 = numpy.attr("float64");

    py::array magnitude_array = py::array::ensure(numpy.attr("asarray")(magnitude, float64));
    if (!magnitude_array || magnitude_array.ndim() != 1) {
        throw py::value_error("radii must be a one dimensional quantity array.");
    }

    py::array_t<double, py::array::c_style | py::array::forcecast> radii_double(magnitude_array);
    auto buf = radii_double.unchecked<1>();

    std::vector<double> out(static_cast<std::size_t>(buf.shape(0)));
    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
        out[static_cast<std::size_t>(i)] = buf(i);
    }
    return out;
}

std::vector<double> array_like_1d_to_double_vector(py::object values) {
    validated_no_unit(values, "number_fractions");
    py::module_ numpy = py::module_::import("numpy");
    py::object float64 = numpy.attr("float64");

    py::array array = py::array::ensure(numpy.attr("asarray")(values, float64));
    if (!array || array.ndim() != 1) {
        throw py::value_error("number_fractions must be a one dimensional array.");
    }

    py::array_t<double, py::array::c_style | py::array::forcecast> as_double(array);
    auto buf = as_double.unchecked<1>();

    std::vector<double> out(static_cast<std::size_t>(buf.shape(0)));
    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
        out[static_cast<std::size_t>(i)] = buf(i);
    }
    return out;
}