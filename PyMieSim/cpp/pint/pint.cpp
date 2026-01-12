#include "pint.h"

namespace py = pybind11;

UnitRegistrySingleton& UnitRegistrySingleton::instance() {
    static UnitRegistrySingleton singleton_instance;
    return singleton_instance;
}

py::object get_shared_ureg() {
    py::gil_scoped_acquire gil;

    py::module_ interface_module = py::module_::import("PackLab.binary.interface_pint");
    if (!py::hasattr(interface_module, "get_ureg")) {
        throw std::runtime_error("PackLab.binary.interface_pint.get_ureg not found.");
    }
    return interface_module.attr("get_ureg")();
}


double to_meters_strict(py::handle value) {
    if (py::isinstance<py::float_>(value) || py::isinstance<py::int_>(value)) {
        throw py::type_error("length must be a pint.Quantity (numbers without units are not accepted)");
    }

    if (!py::hasattr(value, "to")) {
        throw py::type_error("length must be a pint.Quantity (must provide a .to method)");
    }

    py::object q = py::reinterpret_borrow<py::object>(value);
    py::object q_m = q.attr("to")("meter");

    if (!py::hasattr(q_m, "magnitude")) {
        throw py::type_error("length.to('meter') did not return an object with a magnitude attribute");
    }

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
        if (py::isinstance<py::int_>(item) || py::isinstance<py::float_>(item)) {
            throw py::type_error("expected pint.Quantity elements (numbers without units are not accepted)");
        }
        if (!py::hasattr(item, "to")) {
            throw py::type_error("expected pint.Quantity elements (missing .to method)");
        }

        py::object quantity = py::reinterpret_borrow<py::object>(item);
        py::object converted = quantity.attr("to")(unit);

        if (!py::hasattr(converted, "magnitude")) {
            throw py::type_error("quantity.to(unit) did not return an object with a magnitude attribute");
        }

        magnitudes.push_back(py::cast<double>(converted.attr("magnitude")));
    }

    return magnitudes;
}


double quantity_scalar_to_meters(py::object quantity) {
    py::object ureg = get_shared_ureg();
    py::object meter = ureg.attr("meter");
    py::object converted = quantity.attr("to")(meter);
    return py::float_(converted.attr("magnitude"));
}

std::vector<double> quantity_1d_to_meters_vector(py::object quantity) {
    py::object ureg = get_shared_ureg();
    py::object meter = ureg.attr("meter");

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