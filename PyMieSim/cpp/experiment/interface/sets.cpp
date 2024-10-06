#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support

#include "experiment/includes/sets.cpp"

namespace py = pybind11;
typedef std::complex<double> complex128;


PYBIND11_MODULE(SetsInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

// Binding for SPHERE::Set
    py::class_<SPHERE::Set>(module, "CppSphereSet")
        .def(py::init<std::vector<double>, std::vector<complex128>, std::vector<double>>(),
            py::arg("diameter"),
            py::arg("index"),
            py::arg("medium_index"),
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index.")

        .def(py::init<std::vector<double>, std::vector<std::vector<complex128>>, std::vector<double>>(),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_index"),
            "Initializes a set of spheres with given diameters, material indices (for each wavelength), and medium refractive index.")

        .def(py::init<std::vector<double>, std::vector<complex128>, std::vector<std::vector<double>>>(),
            py::arg("diameter"),
            py::arg("index"),
            py::arg("medium_material"),
            "Initializes a set of spheres with given diameters, material indices (for each wavelength), and medium material.")

        .def(py::init<std::vector<double>, std::vector<std::vector<complex128>>, std::vector<std::vector<double>>>(),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_material"),
            "Initializes a set of spheres with given diameters, material indices (for each wavelength), and medium material.");

// Binding for CYLINDER::Set
    py::class_<CYLINDER::Set>(module, "CppCylinderSet")
        .def(py::init<std::vector<double>, std::vector<complex128>, std::vector<double>>(),
            py::arg("diameter"),
            py::arg("index"),
            py::arg("medium_index"),
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index.")

        .def(py::init<std::vector<double>, std::vector<std::vector<complex128>>, std::vector<double>>(),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_index"),
            "Initializes a set of spheres with given diameters, material indices (for each wavelength), and medium refractive index.")

        .def(py::init<std::vector<double>, std::vector<complex128>, std::vector<std::vector<double>>>(),
            py::arg("diameter"),
            py::arg("index"),
            py::arg("medium_material"),
            "Initializes a set of spheres with given diameters, material indices (for each wavelength), and medium material.")

        .def(py::init<std::vector<double>, std::vector<std::vector<complex128>>, std::vector<std::vector<double>>>(),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_material"),
            "Initializes a set of spheres with given diameters, material indices (for each wavelength), and medium material.");

// Binding for CORESHELL::Set
    py::class_<CORESHELL::Set>(module, "CppCoreShellSet")
        .def(py::init<std::vector<double>, std::vector<double>, std::vector<complex128>, std::vector<complex128>, std::vector<double>>(),
            py::arg("core_diameter"),
            py::arg("shell_width"),
            py::arg("core_index"),
            py::arg("shell_index"),
            py::arg("medium_index"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")

        .def(pybind11::init<std::vector<double>&, std::vector<double>&, std::vector<complex128>&, std::vector<std::vector<complex128>>&, std::vector<double>&>(),
            pybind11::arg("core_diameter"),
            pybind11::arg("shell_width"),
            pybind11::arg("core_index"),
            pybind11::arg("shell_material"),
            pybind11::arg("medium_index"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")

        .def(pybind11::init<std::vector<double>&, std::vector<double>&, std::vector<std::vector<complex128>>&, std::vector<complex128>&, std::vector<double>&>(),
            pybind11::arg("core_diameter"),
            pybind11::arg("shell_width"),
            pybind11::arg("core_material"),
            pybind11::arg("shell_index"),
            pybind11::arg("medium_index"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")

        .def(pybind11::init<std::vector<double>&, std::vector<double>&, std::vector<std::vector<complex128>>&, std::vector<std::vector<complex128>>&, std::vector<double>&>(),
            pybind11::arg("core_diameter"),
            pybind11::arg("shell_width"),
            pybind11::arg("core_material"),
            pybind11::arg("shell_material"),
            pybind11::arg("medium_index"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")

        .def(py::init<std::vector<double>, std::vector<double>, std::vector<complex128>, std::vector<complex128>, std::vector<std::vector<double>>>(),
            py::arg("core_diameter"),
            py::arg("shell_width"),
            py::arg("core_index"),
            py::arg("shell_index"),
            py::arg("medium_material"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")

        .def(pybind11::init<std::vector<double>&, std::vector<double>&, std::vector<complex128>&, std::vector<std::vector<complex128>>&, std::vector<std::vector<double>>&>(),
            pybind11::arg("core_diameter"),
            pybind11::arg("shell_width"),
            pybind11::arg("core_index"),
            pybind11::arg("shell_material"),
            pybind11::arg("medium_material"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")

        .def(pybind11::init<std::vector<double>&, std::vector<double>&, std::vector<std::vector<complex128>>&, std::vector<complex128>&, std::vector<std::vector<double>>&>(),
            pybind11::arg("core_diameter"),
            pybind11::arg("shell_width"),
            pybind11::arg("core_material"),
            pybind11::arg("shell_index"),
            pybind11::arg("medium_material"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")

        .def(pybind11::init<std::vector<double>&, std::vector<double>&, std::vector<std::vector<complex128>>&, std::vector<std::vector<complex128>>&, std::vector<std::vector<double>>&>(),
            pybind11::arg("core_diameter"),
            pybind11::arg("shell_width"),
            pybind11::arg("core_material"),
            pybind11::arg("shell_material"),
            pybind11::arg("medium_material"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.");

// Binding for SOURCE::Set
    py::class_<SOURCE::Set>(module, "CppSourceSet")
        .def(py::init<std::vector<double>, std::vector<std::vector<complex128>>, std::vector<double>, std::vector<double>>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            "Initializes a gaussian source set with specific wavelengths, Jones vectors, and amplitudes.")
        .def(py::init<std::vector<double>, std::vector<std::vector<complex128>>, std::vector<double>>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitudes"),
            "Initializes a planewave source set with specific wavelengths, Jones vectors, and amplitudes.");

// Binding for DETECTOR::Set
    py::class_<DETECTOR::Set>(module, "CppDetectorSet")
        .def(py::init<std::vector<std::string>, std::vector<unsigned>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, bool, bool>(),
             py::arg("mode_number"),
             py::arg("sampling"),
             py::arg("NA"),
             py::arg("phi_offset"),
             py::arg("gamma_offset"),
             py::arg("polarization_filter"),
             py::arg("rotation"),
             py::arg("coherent"),
             py::arg("mean_coupling"),
             "Initializes a detector set with scalar fields, numerical aperture, offsets, filters, angle, coherence, and coupling type.");
}

// -
