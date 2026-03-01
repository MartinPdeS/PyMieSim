#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "polarization.h"
#include <pint/pint.h>

namespace py = pybind11;


void register_polarization(py::module_& module)
{
    py::object ureg = get_shared_ureg();

    py::class_<PolarizationState, std::shared_ptr<PolarizationState>>(module, "PolarizationState",
        R"pbdoc(
            Polarization state binding for PyMieSim.
            This class represents the polarization state of the source in PyMieSim.
            It can be initialized either with a Jones vector or with an angle in radians.
        )pbdoc"
        )
        .def(
            py::init(
                [](const std::vector<complex128>& jones_vector) {
                    return std::make_shared<PolarizationState>(jones_vector);
                }
            ),
            py::arg("jones_vector"),
            R"pbdoc(
                Initializes the polarization state using a Jones vector.
                The Jones vector should be a list or array of two complex numbers representing the x and y components of the electric field, respectively.
                Example:
                polarization = PolarizationState(jones_vector=[1+0j, 0+1j])
             )pbdoc"
        )
        .def(py::init<>(
            [](const py::object &angle) {
                double angle_radian = angle.attr("to")("radian").attr("magnitude").cast<double>();
                return std::make_shared<PolarizationState>(angle_radian);
            }),
            py::arg("angle"),
            R"pbdoc(
                Initializes the polarization state using an angle in radians.
                The angle represents the orientation of the polarization, where 0 radians corresponds to linear polarization along the x-axis, and π/4 radians corresponds to circular polarization.


                Parameters
                ----------
                angle : pint.Quantity
                        The angle in radians representing the orientation of the polarization. It can be provided as a pint
             )pbdoc"
        )
        ;

        py::class_<RightCircular, PolarizationState, std::shared_ptr<RightCircular>>(module, "RightCircular",
            R"pbdoc(
                Right circular polarization state.
                This class represents a right circular polarization state, where the electric field rotates in a right-handed manner as it propagates.
                It can be initialized with an optional angle in radians to specify the orientation of the polarization.
            )pbdoc"
        )
        .def(py::init<>())
        ;

        py::class_<LeftCircular, PolarizationState, std::shared_ptr<LeftCircular>>(module, "LeftCircular",
            R"pbdoc(
                Left circular polarization state.
                This class represents a left circular polarization state, where the electric field rotates in a left-handed manner as it propagates.
                It can be initialized with an optional angle in radians to specify the orientation of the polarization.
            )pbdoc"
        )
        .def(py::init<>())
        ;
    }