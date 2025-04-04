#include <pybind11/pybind11.h>

#include "coreshell/coreshell.h"
#include "sphere/sphere.h"
#include "cylinder/cylinder.h"
#include "detector/detector.h"

PYBIND11_MODULE(interface_detector, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<Detector>(module, "BindedDetector")
        .def(py::init<std::string, size_t, double, double, double, double, double, double, bool, bool, double>(),
             py::arg("mode_number"),
             py::arg("sampling"),
             py::arg("NA"),
             py::arg("cache_NA"),
             py::arg("phi_offset"),
             py::arg("gamma_offset"),
             py::arg("polarization_filter"),
             py::arg("rotation"),
             py::arg("coherent"),
             py::arg("mean_coupling"),
             py::arg("medium_refractive_index") = 1.0,
             "Constructs a Detector with given parameters. The `mean_coupling` parameter determines the coupling type (true for point, false for mean).")

        .def("CouplingSphere", &Detector::get_coupling<Sphere>, py::arg("scatterer"), "Calculates the coupling of the detector with a sphere scatterer.")
        .def("CouplingCylinder", &Detector::get_coupling<Cylinder>, py::arg("scatterer"), "Calculates the coupling of the detector with a cylinder scatterer.")
        .def("CouplingCoreShell", &Detector::get_coupling<CoreShell>, py::arg("scatterer"), "Calculates the coupling of the detector with a core-shell scatterer.")
        .def_readwrite("scalar_field", &Detector::scalar_field, "Stores the scalar field values corresponding to the light intensity distribution detected.")
        .def_readwrite("coherent", &Detector::coherent, "Boolean flag indicating whether the detector operates in a coherent detection mode.")
        .def_readonly("NA", &Detector::NA, "Numerical Aperture (NA) of the detector which determines the angular acceptance of light.")
        .def_readonly("sampling", &Detector::sampling, "Samplign of the field.")
        .def_readonly("phi_offset", &Detector::phi_offset, "Offset in the azimuthal angle (phi) used to calibrate the detector orientation.")
        .def_readonly("gamma_offset", &Detector::gamma_offset, "Offset in the polar angle (gamma) used for angular calibration of the detector.")
        .def_readonly("polarization_filter", &Detector::polarization_filter, "Indicates the presence and characteristics of any polarization filter in the detector.")
        .def_readonly("rotation", &Detector::rotation, "The rotation angle of the detector's field of view, typically used in alignment procedures.")
        .def_readonly("mesh", &Detector::fibonacci_mesh, "The Fibonacci mesh used by the detector.")
        .def_readonly("max_angle", &Detector::max_angle, "The Fibonacci mesh max_angle.")
        .def_readonly("min_angle", &Detector::min_angle, "The Fibonacci mesh min_angle [0 if no cache is applied].")
        ;
}

