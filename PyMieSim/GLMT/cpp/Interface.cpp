#include <iostream>
#include <vector>
#include <complex>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h> 
#include <pybind11/numpy.h>
#include "../../includes/SpecialFunc.h"
#include "../../includes/utils.h"
#include "../../includes/BaseClass.h"

namespace py = pybind11;

typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;
typedef py::buffer_info info;

#include "Sphere.cpp"
#include "Cylinder.cpp"


PYBIND11_MODULE(Scatterer, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

    py::class_<_SPHERE>(module, "SPHERE")
    .def(py::init<double, double, double, double, double, double, Cndarray>(),
        py::arg("Index"),
        py::arg("Diameter"),
        py::arg("Wavelength"),
        py::arg("nMedium")      = 1.,
        py::arg("Polarization") = 0.,
        py::arg("E0")           = 1.,
        py::arg("BSC")          = 1. )


     .def("sS1S2",
          &_SPHERE::sS1S2,
          py::arg("Phi"),
          py::arg("Theta"))

     .def("uS1S2",
          &_SPHERE::uS1S2,
          py::arg("Phi"),
          py::arg("Theta"))
    .def("uFields",
         &_SPHERE::uFields,
         py::arg("Phi"),
         py::arg("Theta"),
         py::arg("R") )

    .def("sFields",
         &_SPHERE::sFields,
         py::arg("Phi"),
         py::arg("Theta"),
         py::arg("R") )

    .def("an", &_SPHERE::An, py::arg("MaxOrder")  = 5)

    .def("bn", &_SPHERE::Bn, py::arg("MaxOrder")  = 5)

    .def("cn", &_SPHERE::Cn, py::arg("MaxOrder")  = 5)

    .def("dn", &_SPHERE::Dn, py::arg("MaxOrder")  = 5)

    .def_property_readonly("Efficiencies", &_SPHERE::GetEfficiencies);


  py::class_<_CYLINDER>(module, "CYLINDER")
  .def(py::init<double, double, double, double, double, double, Cndarray>(),
       py::arg("Index"),
       py::arg("Diameter"),
       py::arg("Wavelength"),
       py::arg("nMedium")      = 1.,
       py::arg("Polarization") = 0.,
       py::arg("E0")           = 1.,
       py::arg("BSC")          = 1. )


    .def("sS1S2",
         &_CYLINDER::sS1S2,
         py::arg("Phi"),
         py::arg("Theta"))

    .def("uS1S2",
         &_CYLINDER::uS1S2,
         py::arg("Phi"),
         py::arg("Theta"))

   .def("uFields",
        &_CYLINDER::uFields,
        py::arg("Phi"),
        py::arg("Theta"),
        py::arg("R") )

   .def("sFields",
        &_CYLINDER::sFields,
        py::arg("Phi"),
        py::arg("Theta"),
        py::arg("R") )


   .def("an", &_CYLINDER::An, py::arg("MaxOrder")  = 5)

   .def("bn", &_CYLINDER::Bn, py::arg("MaxOrder")  = 5)

   .def_property_readonly("Efficiencies", &_CYLINDER::GetEfficiencies);




}
















//-
