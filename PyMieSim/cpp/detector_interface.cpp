#include <pybind11/pybind11.h>
#include "special_function.cpp"
#include "detectors.cpp"
#include "cylinder.cpp"
#include "sphere.cpp"
#include "core_shell.cpp"

namespace py = pybind11;
using namespace DETECTOR;

PYBIND11_MODULE(DetectorInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<State>(module, "BindedDetectorState")
        .def_readwrite("scalar_field", &State::scalar_field, "Stores the scalar field values corresponding to the light intensity distribution detected.")
        .def_readwrite("coherent", &State::coherent, "Boolean flag indicating whether the detector operates in a coherent detection mode.")
        .def_readwrite("point_coupling", &State::point_coupling, "Represents the point coupling efficiency of the detector, impacting the detection process.")
        .def_readonly("NA", &State::NA, "Numerical Aperture (NA) of the detector which determines the angular acceptance of light.")
        .def_readonly("sampling", &State::sampling, "Samplign of the field.")
        .def_readonly("phi_offset", &State::phi_offset, "Offset in the azimuthal angle (phi) used to calibrate the detector orientation.")
        .def_readonly("gamma_offset", &State::gamma_offset, "Offset in the polar angle (gamma) used for angular calibration of the detector.")
        .def_readonly("polarization_filter", &State::polarization_filter, "Indicates the presence and characteristics of any polarization filter in the detector.")
        .def_readonly("rotation_angle", &State::rotation_angle, "The rotation angle of the detector's field of view, typically used in alignment procedures.")
        .def_readonly("mesh", &State::fibonacci_mesh, "The Fibonacci mesh used by the detector.");

    // Binding for DETECTOR::Detector class
    py::class_<Detector>(module, "BindedDetector")
        .def(py::init<py::array_t<complex128>, double, double, double, double, double, bool, bool>(),
             py::arg("scalar_field"),
             py::arg("NA"),
             py::arg("phi_offset"),
             py::arg("gamma_offset"),
             py::arg("polarization_filter"),
             py::arg("rotation_angle"),
             py::arg("coherent"),
             py::arg("point_coupling"),
             "Constructs a Detector with given parameters. The `point_coupling` parameter determines the coupling type (true for point, false for mean).")

        .def(py::init<size_t, double, double, double, double, double, bool, bool>(),
             py::arg("sampling"),
             py::arg("NA"),
             py::arg("phi_offset"),
             py::arg("gamma_offset"),
             py::arg("polarization_filter"),
             py::arg("rotation_angle"),
             py::arg("coherent"),
             py::arg("point_coupling"),
             "Constructs a Detector with given parameters. The `point_coupling` parameter determines the coupling type (true for point, false for mean).")

        .def("CouplingSphere", &Detector::get_coupling<SPHERE::Scatterer>, py::arg("scatterer"), "Calculates the coupling of the detector with a sphere scatterer.")
        .def("CouplingCylinder", &Detector::get_coupling<CYLINDER::Scatterer>, py::arg("scatterer"), "Calculates the coupling of the detector with a cylinder scatterer.")
        .def("CouplingCoreShell", &Detector::get_coupling<CORESHELL::Scatterer>, py::arg("scatterer"), "Calculates the coupling of the detector with a core-shell scatterer.")
        .def_readonly("state", &Detector::state, "The state of the detector.");
}
