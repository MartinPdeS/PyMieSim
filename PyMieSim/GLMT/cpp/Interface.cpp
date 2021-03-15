#include <iostream>
#include <vector>
#include <complex>
#include <tuple>
#include <boost/math/special_functions/legendre.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;



//#include "utils.cpp"
#include "Special.cpp"
//#include "Math.cpp"

#include "Sphere.cpp"
#include "Cylinder.cpp"





PYBIND11_MODULE(Cylinder, module) {
    module.doc() = "LGeneralized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

      py::module_ Structured = module.def_submodule("Structured", "Structured fields data");

      py::module_ Unstructured = module.def_submodule("Unstructured", "Unstructured fields data");

    Structured.def("S1S2",
               &Cylinder::S1S2Structured,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("Polarization"),
               py::arg("BSC"),
               py::arg("MaxOrder"),
               "Return S1");

   Structured.def("S1S2Unpolarized",
              &Cylinder::S1S2StructuredUnpolarized,
              py::arg("Index"),
              py::arg("Diameter"),
              py::arg("Wavelength"),
              py::arg("nMedium"),
              py::arg("Phi"),
              py::arg("Theta"),
              py::arg("BSC"),
              py::arg("MaxOrder"),
              "Return S1");


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
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");


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
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");

     Unstructured.def("S1S2",
                &Cylinder::S1S2Unstructured,
                py::arg("Index"),
                py::arg("Diameter"),
                py::arg("Wavelength"),
                py::arg("nMedium"),
                py::arg("Phi"),
                py::arg("Theta"),
                py::arg("Polarization"),
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");

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
               py::arg("BSC"),
               py::arg("MaxOrder"),
               "Return S1");


      Unstructured.def("S1S2Unpolarized",
                 &Cylinder::S1S2UnstructuredUnpolarized,
                 py::arg("Index"),
                 py::arg("Diameter"),
                 py::arg("Wavelength"),
                 py::arg("nMedium"),
                 py::arg("Phi"),
                 py::arg("Theta"),
                 py::arg("BSC"),
                 py::arg("MaxOrder"),
                 "Return S1");


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
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");

}





PYBIND11_MODULE(Sphere, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

    py::module_ Structured = module.def_submodule("Structured", "Structured fields data");

    py::module_ Unstructured = module.def_submodule("Unstructured", "Unstructured fields data");

    Structured.def("S1S2",
               &Sphere::S1S2Structured,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("Polarization"),
               py::arg("BSC"),
               py::arg("MaxOrder"),
               "Return S1");

   Structured.def("S1S2Unpolarized",
              &Sphere::S1S2StructuredUnpolarized,
              py::arg("Index"),
              py::arg("Diameter"),
              py::arg("Wavelength"),
              py::arg("nMedium"),
              py::arg("Phi"),
              py::arg("Theta"),
              py::arg("BSC"),
              py::arg("MaxOrder"),
              "Return S1");


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
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");


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
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");

     Unstructured.def("S1S2",
                &Sphere::S1S2Unstructured,
                py::arg("Index"),
                py::arg("Diameter"),
                py::arg("Wavelength"),
                py::arg("nMedium"),
                py::arg("Phi"),
                py::arg("Theta"),
                py::arg("Polarization"),
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");

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
               py::arg("BSC"),
               py::arg("MaxOrder"),
               "Return S1");


      Unstructured.def("S1S2Unpolarized",
                 &Sphere::S1S2UnstructuredUnpolarized,
                 py::arg("Index"),
                 py::arg("Diameter"),
                 py::arg("Wavelength"),
                 py::arg("nMedium"),
                 py::arg("Phi"),
                 py::arg("Theta"),
                 py::arg("BSC"),
                 py::arg("MaxOrder"),
                 "Return S1");


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
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");



   module.def("an",
              &Sphere::an,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return an");

   module.def("bn",
              &Sphere::bn,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return bn");

   module.def("cn",
              &Sphere::cn,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return cn");

   module.def("dn",
              &Sphere::dn,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return dn");

}

















//-
