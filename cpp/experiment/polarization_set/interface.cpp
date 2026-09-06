#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <complex> // For std::complex support

#include <sstream>

#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include "./polarization_set.h"
#include <utils/casting.h>

using complex128 = std::complex<double>;
namespace py = pybind11;

PYBIND11_MODULE(polarization_set, module) {
    py::object ureg = get_shared_ureg();

    auto format_polarization_set_repr = [](const PolarizationSet& self) {
        std::ostringstream stream;
        stream << "<PolarizationSet states=" << self.number_of_states();

        if (!self.angles.empty()) {
            stream << ", angles=" << self.angles.size();
        }

        stream << ">";
        return stream.str();
    };

    py::class_<PolarizationSet, std::shared_ptr<PolarizationSet>>(
        module,
        "PolarizationSet",
        R"pdoc(
            Ordered collection of polarization states.

            A ``PolarizationSet`` can be created either from polarization angles
            or from Jones vectors. It is used throughout the experiment API to
            describe one or more analyzer or source polarization states.
        )pdoc"
    )
        .def(
            py::init<>(
                [](const py::object& angles) {
                    std::vector<double> angles_value = \
                        Casting::cast_py_to_vector<double>(angles, "radian");
                    return std::make_shared<PolarizationSet>(angles_value);
                }
            ),
            py::arg("angles"),
            R"pdoc(
                Initialize a polarization set from one or more angles.

                Parameters
                ----------
                angles : Angle or array-like of Angle
                    Polarization angles converted to radians internally.
            )pdoc"
        )
        .def(
            py::init<>(
                [](const py::object& angle) {
                    double angle_radian = angle.attr("to")("radian").attr("magnitude").cast<double>();
                    return std::make_shared<PolarizationSet>(angle_radian);
                }
            ),
            py::arg("angles"),
            R"pdoc(
                Initialize a polarization set from a single angle quantity.

                Parameters
                ----------
                angles : Angle
                    Single polarization angle converted to radians internally.
            )pdoc"
        )
        .def(
            py::init(
                [](const std::vector<std::vector<complex128>>& jones_vector) {
                    return std::make_shared<PolarizationSet>(jones_vector);
                }
            ),
            py::arg("jones_vectors"),
            R"pdoc(
                Initialize a polarization set from Jones vectors.

                Parameters
                ----------
                jones_vectors : list[list[complex]]
                    Sequence of Jones vectors, each containing exactly two
                    complex components.
            )pdoc"
        )
        .def(
            "__len__",
            [](const PolarizationSet& self) { return self.number_of_states(); },
            "Returns the number of polarization states defined by the angles."
        )
        .def(
            "__getitem__",
            [ureg](const PolarizationSet& self, size_t index) {
                if (index >= self.number_of_states()) {
                    throw std::out_of_range("Index out of range for PolarizationSet.");
                }
                return py::float_(self.angles[index]) * ureg.attr("radian");
            },
            "Returns the polarization angle at the specified index with appropriate units."
        )
        .def(
            "__iter__",
            [](const PolarizationSet& self) {
                return py::make_iterator(self.angles.begin(), self.angles.end());
            },
            py::keep_alive<0, 1>(), // Keep the PolarizationSet alive while the iterator is used
            "Returns an iterator over the polarization angles in the set."
        )
        .def(
            "__repr__",
            [format_polarization_set_repr](const PolarizationSet& self) {
                return format_polarization_set_repr(self);
            },
            "Return a compact representation showing the number of polarization states."
        )
    ;
}

// -
