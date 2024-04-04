#include <pybind11/pybind11.h>
#include "detectors.cpp"
#include "cylinder.cpp"
#include "sphere.cpp"
#include "core_shell.cpp"

namespace py = pybind11;
using namespace DETECTOR;

PYBIND11_MODULE(DetectorInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    // Binding for DETECTOR::Detector class
    py::class_<Detector>(module, "BindedDetector")
        .def(py::init<CVector, double, double, double, double, double, bool, bool>(),
             py::arg("scalar_field"),
             py::arg("NA"),
             py::arg("phi_offset"),
             py::arg("gamma_offset"),
             py::arg("polarization_filter"),
             py::arg("rotation_angle"),
             py::arg("coherent"),
             py::arg("point_coupling"),
             "Constructs a Detector with given parameters. The `point_coupling` parameter determines the coupling type (true for point, false for mean).")
        .def("rotate_around_axis", &Detector::rotate_around_axis, "Rotates the detector around a specified axis.")
        .def("CouplingSphere", &Detector::get_coupling<SPHERE::Scatterer>, py::arg("scatterer"), "Calculates the coupling of the detector with a sphere scatterer.")
        .def("CouplingCylinder", &Detector::get_coupling<CYLINDER::Scatterer>, py::arg("scatterer"), "Calculates the coupling of the detector with a cylinder scatterer.")
        .def("CouplingCoreShell", &Detector::get_coupling<CORESHELL::Scatterer>, py::arg("scatterer"), "Calculates the coupling of the detector with a core-shell scatterer.")
        .def_readwrite("mesh", &Detector::fibonacci_mesh, "The Fibonacci mesh used by the detector.");
}
