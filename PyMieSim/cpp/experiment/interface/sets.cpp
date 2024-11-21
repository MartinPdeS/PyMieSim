#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support

#include "experiment/includes/scatterer_properties.cpp"
#include "experiment/includes/sets.cpp"

namespace py = pybind11;
typedef std::complex<double> complex128;


PYBIND11_MODULE(SetsInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

// Binding for SPHERE::Set
    py::class_<SPHERE::Set>(module, "CppSphereSet")
        .def(py::init<const std::vector<double>&, const ScattererProperties&, const MediumProperties&>(),
            py::arg("diameter"),
            py::arg("scatterer_properties"),
            py::arg("medium_properties"),
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index.")
            ;

// Binding for CYLINDER::Set
    py::class_<CYLINDER::Set>(module, "CppCylinderSet")
        .def(py::init<const std::vector<double>&, const ScattererProperties&, const MediumProperties&>(),
            py::arg("diameter"),
            py::arg("scatterer_properties"),
            py::arg("medium_properties"),
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index.")
            ;

// Binding for CORESHELL::Set
    py::class_<CORESHELL::Set>(module, "CppCoreShellSet")
        .def(py::init<const std::vector<double>&, const std::vector<double>&, const ScattererProperties&, const ScattererProperties&, const MediumProperties&>(),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_properties"),
            py::arg("shell_properties"),
            py::arg("medium_properties"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")
            ;

// Binding for SOURCE::Set
    py::class_<SOURCE::Set>(module, "CppSourceSet")
        .def(py::init<const std::vector<double>&, const std::vector<std::vector<complex128>>&, const std::vector<double>&, const std::vector<double>&>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            "Initializes a gaussian source set with specific wavelengths, Jones vectors, and amplitudes.")
        .def(py::init<const std::vector<double>&, const std::vector<std::vector<complex128>>&, const std::vector<double>&>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            "Initializes a planewave source set with specific wavelengths, Jones vectors, and amplitudes.");

// Binding for DETECTOR::Set
    py::class_<DETECTOR::Set>(module, "CppDetectorSet")
        .def(py::init<const std::vector<std::string>&, const std::vector<unsigned>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, bool, bool>(),
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
             "Initializes a detector set with scalar fields, numerical aperture, offsets, filters, angle, coherence, and coupling type.");
}

// -
