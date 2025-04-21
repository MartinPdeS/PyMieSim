#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support

#include "sets/properties.h"
#include "sets/source.h"
#include "sets/scatterer.h"
#include "sets/detector.h"

namespace py = pybind11;
typedef std::complex<double> complex128;


PYBIND11_MODULE(interface_sets, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<ScattererProperties>(module, "CppScattererProperties")
        .def(py::init<std::vector<complex128>>(), py::arg("index_properties"))
        .def(py::init<std::vector<std::vector<complex128>>>(), py::arg("material_properties"))
        ;

    py::class_<MediumProperties>(module, "CppMediumProperties")
        .def(py::init<std::vector<double>>(), py::arg("index_properties"))
        .def(py::init<std::vector<std::vector<double>>>(), py::arg("material_properties"))
        ;

    // Binding for SPHERE::Set
    py::class_<SphereSet>(module, "CppSphereSet")
        .def(py::init<const std::vector<double>&, const ScattererProperties&, const MediumProperties&, bool>(),
            py::arg("diameter"),
            py::arg("scatterer_properties"),
            py::arg("medium_properties"),
            py::arg("is_sequential"),
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index.")
            ;

    // Binding for CYLINDER::Set
    py::class_<CylinderSet>(module, "CppCylinderSet")
        .def(py::init<const std::vector<double>&, const ScattererProperties&, const MediumProperties&, bool>(),
            py::arg("diameter"),
            py::arg("scatterer_properties"),
            py::arg("medium_properties"),
            py::arg("is_sequential"),
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index.")
            ;

    // Binding for CORESHELL::Set
    py::class_<CoreShellSet>(module, "CppCoreShellSet")
        .def(py::init<const std::vector<double>&, const std::vector<double>&, const ScattererProperties&, const ScattererProperties&, const MediumProperties&, bool>(),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_properties"),
            py::arg("shell_properties"),
            py::arg("medium_properties"),
            py::arg("is_sequential"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")
            ;

    // Binding for SOURCE::Set
    py::class_<BaseSourceSet>(module, "CppBaseSourceSet");

    py::class_<GaussianSourceSet, BaseSourceSet>(module, "CppGaussianSourceSet")
        .def(py::init<const std::vector<double>&, const std::vector<std::vector<complex128>>&, const std::vector<double>&, const std::vector<double>&, bool>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            py::arg("is_sequential"),
            "Initializes a gaussian source set with specific wavelengths, Jones vectors, and amplitudes.")
        ;

    py::class_<PlaneWaveSourceSet, BaseSourceSet>(module, "CppPlaneWaveSourceSet")
        .def(py::init<const std::vector<double>&, const std::vector<std::vector<complex128>>&, const std::vector<double>&, bool>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            py::arg("is_sequential"),
            "Initializes a planewave source set with specific wavelengths, Jones vectors, and amplitudes.");

    // Binding for DETECTOR::Set
    py::class_<DetectorSet>(module, "CppDetectorSet")
        .def(py::init<>())
        .def(py::init<const std::vector<std::string>&, const std::vector<unsigned>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, bool, bool, bool>(),
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
             py::arg("is_sequential"),
             "Initializes a detector set with scalar fields, numerical aperture, offsets, filters, angle, coherence, and coupling type.");
}

// -