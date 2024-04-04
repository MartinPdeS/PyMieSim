#include <pybind11/pybind11.h>
#include "sources.cpp"  // Assuming this includes the necessary DETECTOR definitions

namespace py = pybind11;
using DETECTOR::Detector;

PYBIND11_MODULE(DetectorInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<Detector>(module, "BindedDetector")
        .def(py::init<CVector, double, double, double, double, bool, bool>(),
             py::arg("ScalarField"),
             py::arg("NA"),
             py::arg("PhiOffset"),
             py::arg("GammaOffset"),
             py::arg("Filter"),
             py::arg("Coherent"),
             py::arg("PointCoupling"),
             "Constructs a Detector with specified optical properties and detection parameters. "
             "`PointCoupling` indicates whether to use point (true) or mean (false) coupling.")

        .def("CouplingSphere", &Detector::Coupling<SPHERE::Scatterer>, py::arg("Scatterer"), "Calculates coupling with a sphere scatterer.")
        .def("CouplingCylinder", &Detector::Coupling<CYLINDER::Scatterer>, py::arg("Scatterer"), "Calculates coupling with a cylinder scatterer.")
        .def("CouplingCoreShell", &Detector::Coupling<CORESHELL::Scatterer>, py::arg("Scatterer"), "Calculates coupling with a core-shell scatterer.");
}
