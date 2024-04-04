#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support
#include "sets.cpp"

namespace py = pybind11;


PYBIND11_MODULE(Sets, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    // Binding for SPHERE::Set
    py::class_<SPHERE::Set>(module, "CppSphereSet")
        .def(py::init<DVector, std::vector<std::complex<double>>, DVector>(),
             py::arg("diameter"),
             py::arg("index"),
             py::arg("n_medium"),
             "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index.")
        .def(py::init<DVector, std::vector<std::vector<std::complex<double>>>, DVector>(),
             py::arg("diameter"),
             py::arg("material_index"),
             py::arg("n_medium"),
             "Initializes a set of spheres with given diameters, material indices (for each wavelength), and medium refractive index.");

    // Binding for CYLINDER::Set
    py::class_<CYLINDER::Set>(module, "CppCylinderSet")
        .def(py::init<DVector, CVector, DVector>(),
             py::arg("diameter"),
             py::arg("index"),
             py::arg("n_medium"),
             "Initializes a set of cylinders with given diameters, refractive indices, and medium refractive index.")
        .def(py::init<DVector, std::vector<std::vector<std::complex<double>>>, DVector>(),
             py::arg("diameter"),
             py::arg("material_index"),
             py::arg("n_medium"),
             "Initializes a set of cylinders with given diameters, material indices (for each wavelength), and medium refractive index.");

    // Binding for CORESHELL::Set
    py::class_<CORESHELL::Set>(module, "CppCoreShellSet")
        .def(py::init<DVector, DVector, CVector, CVector, DVector>(),
             py::arg("core_diameter"),
             py::arg("shell_width"),
             py::arg("core_index"),
             py::arg("shell_index"),
             py::arg("n_medium"),
             "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")
         .def(pybind11::init<DVector&, DVector&, CVector&, std::vector<CVector>&, DVector&>(),
               pybind11::arg("core_diameter"),
               pybind11::arg("shell_width"),
               pybind11::arg("core_index"),
               pybind11::arg("shell_material_index"),
               pybind11::arg("n_medium"),
               "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")
         .def(pybind11::init<DVector&, DVector&, std::vector<CVector>&, CVector&, DVector&>(),
               pybind11::arg("core_diameter"),
               pybind11::arg("shell_width"),
               pybind11::arg("core_material_index"),
               pybind11::arg("shell_index"),
               pybind11::arg("n_medium"),
               "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.")
         .def(pybind11::init<DVector&, DVector&, std::vector<CVector>&, std::vector<CVector>&, DVector&>(),
               pybind11::arg("core_diameter"),
               pybind11::arg("shell_width"),
               pybind11::arg("core_material_index"),
               pybind11::arg("shell_material_index"),
               pybind11::arg("n_medium"),
               "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index.");

    // Binding for SOURCE::Set
    py::class_<SOURCE::Set>(module, "CppSourceSet")
        .def(py::init<DVector, std::vector<CVector>, DVector>(),
             py::arg("wavelength"),
             py::arg("jones_vector"),
             py::arg("amplitude"),
             "Initializes a source set with specific wavelengths, Jones vectors, and amplitudes.");

    // Binding for DETECTOR::Set
    py::class_<DETECTOR::Set>(module, "CppDetectorSet")
        .def(py::init<std::vector<CVector>, DVector, DVector, DVector, DVector, DVector, bool, bool>(),
             py::arg("scalarfield"),
             py::arg("NA"),
             py::arg("phi_offset"),
             py::arg("gamma_offset"),
             py::arg("polarization_filter"),
             py::arg("rotation_angle"),
             py::arg("coherent"),
             py::arg("point_coupling"),
             "Initializes a detector set with scalar fields, numerical aperture, offsets, filters, angle, coherence, and coupling type.");
}

// -
