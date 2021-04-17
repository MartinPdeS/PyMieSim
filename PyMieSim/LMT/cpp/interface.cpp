#include "../includes/SpecialFunc.h"
#include "../includes/utils.h"
#include "../includes/BaseClass.h"
#include "../includes/BaseFunc.h"
#include "Sphere.cpp"
#include "Cylinder.cpp"
#include <iostream>





PYBIND11_MODULE(Scatterer, module) {
    module.doc() = "Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";


      py::class_<SPHERE>(module, "SPHERE")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("Index"),
           py::arg("Diameter"),
           py::arg("Wavelength"),
           py::arg("nMedium")      = 1.,
           py::arg("Polarization") = 0.,
           py::arg("E0")           = 1. )

       .def("S1S2",
            &SPHERE::S1S2,
            py::arg("Phi") )

      .def("uFields",
           &SPHERE::uFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("sFields",
           &SPHERE::sFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("uS1S2",
           &SPHERE::uS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("sS1S2",
           &SPHERE::sS1S2,
           py::arg("Phi"),
           py::arg("Theta")) 

      .def("an", &SPHERE::An, py::arg("MaxOrder")  = 5)

      .def("bn", &SPHERE::Bn, py::arg("MaxOrder")  = 5)

      .def("cn", &SPHERE::Cn, py::arg("MaxOrder")  = 5)

      .def("dn", &SPHERE::Dn, py::arg("MaxOrder")  = 5)

      .def_property_readonly("Efficiencies", &SPHERE::GetEfficiencies);


      py::class_<CYLINDER>(module, "CYLINDER")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("Index"),
           py::arg("Diameter"),
           py::arg("Wavelength"),
           py::arg("nMedium")      = 1.,
           py::arg("Polarization") = 0.,
           py::arg("E0")           = 1. )

       .def("S1S2",
            &CYLINDER::S1S2,
            py::arg("Phi") )

      .def("uFields",
           &CYLINDER::uFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("sFields",
           &CYLINDER::sFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("uS1S2",
           &CYLINDER::uS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("sS1S2",
           &CYLINDER::sS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("an", &CYLINDER::An, py::arg("MaxOrder")  = 5)

      .def("bn", &CYLINDER::Bn, py::arg("MaxOrder")  = 5)


      .def_property_readonly("Efficiencies", &CYLINDER::GetEfficiencies);
}







// -
