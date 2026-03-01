#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <complex> // For std::complex support

#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include "./polarization_set.h"

using complex128 = std::complex<double>;
namespace py = pybind11;

void register_polarization_set(py::module& module) {
    py::object ureg = get_shared_ureg();

    py::class_<PolarizationSet, std::shared_ptr<PolarizationSet>>(module, "PolarizationSet")
        .def(
            py::init<>(
                [](const py::object& angles) {
                    std::vector<double> angles_value = \
                        cast_scalar_or_array_to_vector_double(angles.attr("to")("radian").attr("magnitude"));
                    return std::make_shared<PolarizationSet>(angles_value);
                }
            ),
            py::arg("angles")
        )
        .def(
            py::init<>(
                [](const py::object& angle) {
                    double angle_radian = angle.attr("to")("radian").attr("magnitude").cast<double>();
                    return std::make_shared<PolarizationSet>(angle_radian);
                }
            ),
            py::arg("angles")
        )
        .def(
            py::init(
                [](const std::vector<std::vector<complex128>>& jones_vector) {
                    return std::make_shared<PolarizationSet>(jones_vector);
                }
            ),
            py::arg("jones_vectors")
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
    ;
}

// -
