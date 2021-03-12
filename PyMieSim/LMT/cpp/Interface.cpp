#include "Math.cpp"
#include "Functions.cpp"
#include "Sphere.cpp"
#include "Cylinder.cpp"
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

typedef std::complex<double> complex128;

typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray ;





PYBIND11_MODULE(Sphere, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) c++ binding module for light scattering from a spherical scatterer";

    py::module_ Structured = module.def_submodule("Structured", "Structured fields data");

    py::module_ Unstructured = module.def_submodule("Unstructured", "Unstructured fields data");

    Unstructured.def("Fields",
               &Sphere::FieldsUnstructured,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"), 
               py::arg("Theta"),
               py::arg("Polarization"),
               py::arg("E0"),
               py::arg("R"),
               "Compute the scattering far-field for a spherical scatterer");


   Structured.def("Fields",
              &Sphere::FieldsStructured,
              py::arg("Index"),
              py::arg("Diameter"),
              py::arg("Wavelength"),
              py::arg("nMedium"),
              py::arg("Phi"),
              py::arg("Theta"),
              py::arg("Polarization"),
              py::arg("E0"),
              py::arg("R"),
              "Compute the scattering far-field for a spherical scatterer");


    Unstructured.def("FieldsUnpolarized",
               &Sphere::FieldsUnstructuredUnpolarized,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("E0"),
               py::arg("R"),
               "Compute the scattering far-field for a spherical scatterer");


   Structured.def("FieldsUnpolarized",
              &Sphere::FieldsStructuredUnpolarized,
              py::arg("Index"),
              py::arg("Diameter"),
              py::arg("Wavelength"),
              py::arg("nMedium"),
              py::arg("Phi"),
              py::arg("Theta"),
              py::arg("E0"),
              py::arg("R"),
              "Compute the scattering far-field for a spherical scatterer");


     module.def("S1S2",
                &Sphere::S1S2,
                py::arg("Index"),
                py::arg("Diameter"),
                py::arg("Wavelength"),
                py::arg("nMedium"),
                py::arg("Phi"),
                "Compute the scattering coefficient S1 & S2");


    module.def("Efficiencies",
               &Sphere::Efficiencies,
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("Index"),
               py::arg("nMedium"),
               "Compute the scattering efficiencies");

}



PYBIND11_MODULE(Cylinder, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) c++ binding module for light scattering from a spherical scatterer";

    py::module_ Structured = module.def_submodule("Structured", "Structured fields data");

    py::module_ Unstructured = module.def_submodule("Unstructured", "Unstructured fields data");

    Unstructured.def("Fields",
               &Cylinder::FieldsUnstructured,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("Polarization"),
               py::arg("E0"),
               py::arg("R"),
               "Compute the scattering far-field for a spherical scatterer");


   Structured.def("Fields",
              &Cylinder::FieldsStructured,
              py::arg("Index"),
              py::arg("Diameter"),
              py::arg("Wavelength"),
              py::arg("nMedium"),
              py::arg("Phi"),
              py::arg("Theta"),
              py::arg("Polarization"),
              py::arg("E0"),
              py::arg("R"),
              "Compute the scattering far-field for a spherical scatterer");


    Unstructured.def("FieldsUnpolarized",
               &Cylinder::FieldsUnstructuredUnpolarized,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("E0"),
               py::arg("R"),
               "Compute the scattering far-field for a spherical scatterer");


   Structured.def("FieldsUnpolarized",
              &Cylinder::FieldsStructuredUnpolarized,
              py::arg("Index"),
              py::arg("Diameter"),
              py::arg("Wavelength"),
              py::arg("nMedium"),
              py::arg("Phi"),
              py::arg("Theta"),
              py::arg("E0"),
              py::arg("R"),
              "Compute the scattering far-field for a spherical scatterer");


     module.def("S1S2",
                &Cylinder::S1S2,
                py::arg("Index"),
                py::arg("Diameter"),
                py::arg("Wavelength"),
                py::arg("nMedium"),
                py::arg("Phi"),
                "Compute the scattering coefficient S1 & S2");


      module.def("Efficiencies",
                 &Cylinder::Efficiencies,
                 py::arg("Diameter"),
                 py::arg("Wavelength"),
                 py::arg("Index"),
                 py::arg("nMedium"),
                 "Compute the scattering efficiencies");

}







// -
