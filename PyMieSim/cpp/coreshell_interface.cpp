#include <pybind11/pybind11.h>
#include "core_shell.cpp"

PYBIND11_MODULE(CoreShellInterface, module)
{
     module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

     pybind11::class_<CORESHELL::State>(module, "CppCoreShellState");

     pybind11::class_<CORESHELL::Scatterer>(module, "CORESHELL")
          .def(
               pybind11::init<double &, double &, double &, double &, complex128 &, complex128 &, double &, CVector &>(),
               pybind11::arg("wavelength"),
               pybind11::arg("amplitude"),
               pybind11::arg("core_diameter"),
               pybind11::arg("shell_width"),
               pybind11::arg("core_index"),
               pybind11::arg("shell_index"),
               pybind11::arg("n_medium"),
               pybind11::arg("jones_vector")
          )

          .def(
               "get_s1s2",
               &CORESHELL::Scatterer::get_s1s2_py,
               pybind11::arg("phi")
          )

          .def(
               "get_fields",
               &CORESHELL::Scatterer::get_unstructured_fields_py,
               pybind11::arg("phi"),
               pybind11::arg("theta"),
               pybind11::arg("r")
          )

          .def(
               "get_full_fields",
               &CORESHELL::Scatterer::get_full_structured_fields_py,
               pybind11::arg("sampling"),
               pybind11::arg("r")
          )

     .def("an", pybind11::overload_cast<>(&CORESHELL::Scatterer::get_an_py))
     .def("bn", pybind11::overload_cast<>(&CORESHELL::Scatterer::get_bn_py))

     .def_property_readonly("Qsca", &CORESHELL::Scatterer::get_Qsca)
     .def_property_readonly("Qext", &CORESHELL::Scatterer::get_Qext)
     .def_property_readonly("Qabs", &CORESHELL::Scatterer::get_Qabs)
     .def_property_readonly("Qback", &CORESHELL::Scatterer::get_Qback)
     .def_property_readonly("Qpr", &CORESHELL::Scatterer::get_Qpr)

     .def_property_readonly("Csca", &CORESHELL::Scatterer::get_Csca)
     .def_property_readonly("Cext", &CORESHELL::Scatterer::get_Cext)
     .def_property_readonly("Cabs", &CORESHELL::Scatterer::get_Cabs)
     .def_property_readonly("Cback", &CORESHELL::Scatterer::get_Cback)
     .def_property_readonly("Cpr", &CORESHELL::Scatterer::get_Cpr)

     .def_property_readonly("g", &CORESHELL::Scatterer::get_g)

     .def_readwrite("state", &CORESHELL::Scatterer::state)
     .def_readwrite("area", &CORESHELL::Scatterer::area)
     .def_readwrite("size_parameter", &CORESHELL::Scatterer::size_parameter)
     ;
}







// -
