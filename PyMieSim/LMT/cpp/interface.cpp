//#include "../F90/include/complex_bessel.h"
#include "../includes/SpecialFunc.h"
#include "../includes/utils.h"
#include "../includes/BaseFunc.h"
#include "../includes/BaseClass.h"
#include "Sphere.cpp"
#include "ShellSphere1.cpp"
#include "Cylinder.cpp"
#include <iostream>





PYBIND11_MODULE(Scatterer, module) {
    module.doc() = "Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";


      py::class_<SPHERE>(module, "SPHERE")
      .def(py::init<complex128, double, double, double, double, double>(),
           py::arg("Index"),
           py::arg("Diameter"),
           py::arg("Wavelength"),
           py::arg("nMedium")      = 1.,
           py::arg("Polarization") = 0.,
           py::arg("E0")           = 1. )

       .def("S1S2",
            &BASE::S1S2,
            py::arg("Phi") )

      .def("uFields",
           &BASE::uFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("sFields",
           &BASE::sFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("uS1S2",
           &BASE::uS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("sS1S2",
           &BASE::sS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("an", &SPHERE::An, py::arg("MaxOrder")  = 5)

      .def("bn", &SPHERE::Bn, py::arg("MaxOrder")  = 5)

      .def("cn", &SPHERE::Cn, py::arg("MaxOrder")  = 5)

      .def("dn", &SPHERE::Dn, py::arg("MaxOrder")  = 5)

      .def_property_readonly("Efficiencies", &BASE::GetEfficiencies);




      py::class_<SHELLSPHERE1>(module, "SHELLSPHERE1")
      .def(py::init<complex128,complex128, double, double, double, double, double, double>(),
           py::arg("ShellIndex"),
           py::arg("CoreIndex"),
           py::arg("ShellDiameter"),
           py::arg("CoreDiameter"),
           py::arg("Wavelength"),
           py::arg("nMedium")      = 1.,
           py::arg("Polarization") = 0.,
           py::arg("E0")           = 1. )

       .def("S1S2",
            &BASE::S1S2,
            py::arg("Phi") )

      .def("uFields",
           &BASE::uFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("sFields",
           &BASE::sFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("uS1S2",
           &BASE::uS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("sS1S2",
           &BASE::sS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("an", &SHELLSPHERE1::An, py::arg("MaxOrder")  = 5)

      .def("bn", &SHELLSPHERE1::Bn, py::arg("MaxOrder")  = 5)

      .def_property_readonly("Efficiencies", &SHELLSPHERE1::GetEfficiencies);






      py::class_<CYLINDER>(module, "CYLINDER")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("Index"),
           py::arg("Diameter"),
           py::arg("Wavelength"),
           py::arg("nMedium")      = 1.,
           py::arg("Polarization") = 0.,
           py::arg("E0")           = 1. )

       .def("S1S2",
            &BASE::S1S2,
            py::arg("Phi") )

      .def("uFields",
           &BASE::uFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("sFields",
           &BASE::sFields,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R") )

      .def("uS1S2",
           &BASE::uS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("sS1S2",
           &BASE::sS1S2,
           py::arg("Phi"),
           py::arg("Theta"))

      .def("an", &CYLINDER::An, py::arg("MaxOrder")  = 5)

      .def("bn", &CYLINDER::Bn, py::arg("MaxOrder")  = 5)

      .def_property_readonly("Efficiencies", &BASE::GetEfficiencies);
}







// -
