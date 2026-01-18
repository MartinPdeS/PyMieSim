#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <complex>
#include <vector>
#include <array>
#include <stdexcept>
#include <string>

#include "polarization.h"
#include <pint/pint.h>

namespace py = pybind11;


using Complex = std::complex<double>;
using Row = std::array<Complex, 2>;

static py::array_t<Complex> element_to_numpy(const JonesVector::Element& element)
{
    const py::ssize_t n_rows = static_cast<py::ssize_t>(element.size());
    py::array_t<Complex> out({n_rows, static_cast<py::ssize_t>(2)});

    auto out_view = out.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < n_rows; ++i) {
        out_view(i, 0) = element[static_cast<std::size_t>(i)][0];
        out_view(i, 1) = element[static_cast<std::size_t>(i)][1];
    }

    return out;
}

static JonesVector::Element numpy_to_element(py::handle obj)
{
    py::array arr_any = py::array::ensure(obj);
    if (!arr_any) {
        throw std::invalid_argument("element must be array like.");
    }

    py::array_t<Complex, py::array::c_style | py::array::forcecast> arr(arr_any);

    if (arr.ndim() == 1) {
        if (arr.shape(0) != 2) {
            throw std::invalid_argument("If element is 1D, it must have length 2.");
        }

        auto r = arr.unchecked<1>();
        JonesVector::Element element;
        element.push_back(Row{r(0), r(1)});
        return element;
    }

    if (arr.ndim() != 2) {
        throw std::invalid_argument("element must be a 1D (length 2) or 2D (N, 2) array.");
    }

    if (arr.shape(1) != 2) {
        throw std::invalid_argument("element must have shape (N, 2).");
    }

    const std::size_t n_rows = static_cast<std::size_t>(arr.shape(0));
    auto r = arr.unchecked<2>();

    JonesVector::Element element;
    element.reserve(n_rows);

    for (std::size_t i = 0; i < n_rows; ++i) {
        element.push_back(Row{r(static_cast<py::ssize_t>(i), 0), r(static_cast<py::ssize_t>(i), 1)});
    }

    return element;
}

static std::vector<double> quantity_to_magnitudes(py::object quantity)
{
    py::object np = py::module_::import("numpy");
    py::object atleast_1d = np.attr("atleast_1d");

    py::object mag = quantity.attr("magnitude");
    py::object mag_1d = atleast_1d(mag);

    py::array mag_any = py::array::ensure(mag_1d);
    if (!mag_any) {
        throw std::invalid_argument("Angle magnitude must be array like.");
    }

    py::array_t<double, py::array::c_style | py::array::forcecast> mag_arr(mag_any);

    std::vector<double> out;
    out.resize(static_cast<std::size_t>(mag_arr.shape(0)));

    auto r = mag_arr.unchecked<1>();
    for (py::ssize_t i = 0; i < mag_arr.shape(0); ++i) {
        out[static_cast<std::size_t>(i)] = r(i);
    }

    return out;
}

static AngleUnit detect_angle_unit_from_quantity(py::object quantity)
{
    py::object ureg = get_shared_ureg();
    py::object units = quantity.attr("units");

    if (units.equal(ureg.attr("degree"))) {
        return AngleUnit::Degree;
    }
    if (units.equal(ureg.attr("radian"))) {
        return AngleUnit::Radian;
    }

    // Fallback heuristic: try converting to degree; if it works, treat as degree based input unit.
    try {
        (void)quantity.attr("to")(ureg.attr("degree"));
        return AngleUnit::Degree;
    } catch (...) {
    }

    // Try converting to radian
    try {
        (void)quantity.attr("to")(ureg.attr("radian"));
        return AngleUnit::Radian;
    } catch (...) {
    }

    throw std::invalid_argument("Angle must be convertible to degree or radian.");
}

static py::object build_angle_quantity(const std::vector<double>& magnitudes, AngleUnit unit)
{
    py::object ureg = get_shared_ureg();
    py::object np = py::module_::import("numpy");

    py::array_t<double> mag_arr(static_cast<py::ssize_t>(magnitudes.size()));
    auto w = mag_arr.mutable_unchecked<1>();
    for (py::ssize_t i = 0; i < mag_arr.shape(0); ++i) {
        w(i) = magnitudes[static_cast<std::size_t>(i)];
    }

    py::object unit_obj = (unit == AngleUnit::Degree) ? ureg.attr("degree") : ureg.attr("radian");
    return mag_arr * unit_obj;
}

static py::object polarization_add(const JonesVector& self, const JonesVector& other)
{
    if (self.has_angle() && other.has_angle()) {
        if (self.angle_unit() != other.angle_unit()) {
            throw std::invalid_argument("Cannot add polarizations with different angle units.");
        }

        std::vector<double> merged = self.angles();
        const auto& b = other.angles();
        merged.insert(merged.end(), b.begin(), b.end());

        Linear result(std::move(merged), self.angle_unit());
        return py::cast(std::move(result));
    }

    JonesVector result = self + other;
    return py::cast(std::move(result));
}












void register_polarization(py::module_& module)
{
    py::object ureg = get_shared_ureg();

    py::enum_<AngleUnit>(module, "AngleUnit")
        .value("Degree", AngleUnit::Degree)
        .value("Radian", AngleUnit::Radian);

    py::class_<BasePolarization>(module, "BasePolarization")
        .def(py::init<>());

    py::class_<JonesVector, BasePolarization>(module, "JonesVector")
        .def(
            py::init(
                [ureg](py::object element) {
                    auto rows = numpy_to_element(element);
                    return JonesVector(std::move(rows));
                }
        ),
        py::arg("element"),
        R"pbdoc(
            Initialize a JonesVector polarization from its elements.

            Parameters
            ----------
            element : array-like, shape (2,) or (N, 2)
                The Jones vector elements as a 1D array of length 2 (single vector)
                or a 2D array of shape (N, 2) (multiple vectors).
        )pbdoc"
        )
        .def_property_readonly(
            "element",
            [ureg](const JonesVector& self) {
                return element_to_numpy(self.elements());
            },
            "The Jones vector elements as a NumPy array."
        )
        .def(
            "__iter__",
            [ureg](const JonesVector& self) {
                if (self.has_angle())
                    return py::iter(build_angle_quantity(self.angles(), self.angle_unit()));

                return py::iter(element_to_numpy(self.elements()));
            },
            "Iterator over the polarization data."
        )
        .def(
            "__add__",
            [ureg](const JonesVector& self, const JonesVector& other) {
                return polarization_add(self, other);
            },
            "Add two JonesVector polarizations together."
        )
        .def(
            "__repr__",
            [ureg](const JonesVector& self) {
                return py::str(self.to_string());
            },
            R"pbdoc(
                Official string representation of the polarization.
            )pbdoc"
        )

        .def(
            "__str__",
            [ureg](const JonesVector& self) {
                return py::str(self.to_string());
            },
            "String representation of the polarization."
        );

    py::class_<RightCircular, JonesVector>(module, "RightCircular")
        .def(py::init<>());

    py::class_<LeftCircular, JonesVector>(module, "LeftCircular")
        .def(py::init<>());

    py::class_<Linear, JonesVector>(module, "Linear")
        .def(py::init([ureg](py::object angle_quantity) {
                const AngleUnit unit = detect_angle_unit_from_quantity(angle_quantity);
                const std::vector<double> magnitudes = quantity_to_magnitudes(angle_quantity);
                return Linear(magnitudes, unit);
            }
        ),
        py::arg("element"))
        .def_property_readonly(
            "angle",
            [ureg](const JonesVector& self) -> py::object {
                if (!self.has_angle())
                    throw py::attribute_error("This polarization has no attribute 'angle'.");

                return build_angle_quantity(self.angles(), self.angle_unit());
            },
            "The polarization angles as a quantity array."
        )
        .def_property_readonly(
            "element",
            [ureg](const JonesVector& self) {
                return element_to_numpy(self.elements());
            },
            "The Jones vector elements as a NumPy array."
        )
        .def(
            "__iter__",
            [ureg](const JonesVector& self) {
                if (self.has_angle())
                    return py::iter(build_angle_quantity(self.angles(), self.angle_unit()));

                return py::iter(element_to_numpy(self.elements()));
            },
            "Iterator over the polarization data."
        )
        .def(
            "__add__",
            [ureg](const JonesVector& self, const JonesVector& other) {
                return polarization_add(self, other);
            },
            "Add two JonesVector polarizations together."
        )
        .def(
            "__repr__",
            [ureg](const JonesVector& self) {
                return py::str(self.to_string());
            },
            "Official string representation of the polarization."
        )
        .def(
            "__str__",
            [ureg](const JonesVector& self) {
                return py::str(self.to_string());
            },
            "String representation of the polarization."
        );
}
