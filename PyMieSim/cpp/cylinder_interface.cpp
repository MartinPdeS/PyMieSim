#include <pybind11/pybind11.h>
#include "cylinder.cpp"

PYBIND11_MODULE(CylinderInterface, module)
{
     module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";


     pybind11::class_<CYLINDER::State>(module, "CppCylinderState");


     pybind11::class_<CYLINDER::Scatterer>(module, "CYLINDER")
          .def(
               pybind11::init<double&, double&, double&, complex128&, double&, CVector&>(),
               pybind11::arg("wavelength"),
               pybind11::arg("amplitude"),
               pybind11::arg("diameter"),
               pybind11::arg("index"),
               pybind11::arg("n_medium"),
               pybind11::arg("jones_vector")
          )

          .def(
               "get_s1s2",
               &CYLINDER::Scatterer::get_s1s2_py,
               pybind11::arg("phi")
          )

          .def(
               "get_fields",
               &CYLINDER::Scatterer::get_unstructured_fields_py,
               pybind11::arg("phi"),
               pybind11::arg("theta"),
               pybind11::arg("r")
          )

          .def(
               "get_full_fields",
               &CYLINDER::Scatterer::get_full_structured_fields_py,
               pybind11::arg("sampling"),
               pybind11::arg("r")
          )

      .def("a1n", pybind11::overload_cast<>(&CYLINDER::Scatterer::get_a1n_py))
      .def("b1n", pybind11::overload_cast<>(&CYLINDER::Scatterer::get_b1n_py))
      .def("a2n", pybind11::overload_cast<>(&CYLINDER::Scatterer::get_a2n_py))
      .def("b2n", pybind11::overload_cast<>(&CYLINDER::Scatterer::get_b2n_py))

      .def_property_readonly("Qsca",  &CYLINDER::Scatterer::get_Qsca)
      .def_property_readonly("Qext",  &CYLINDER::Scatterer::get_Qext)
      .def_property_readonly("Qabs",  &CYLINDER::Scatterer::get_Qabs)
      .def_property_readonly("Qback", &CYLINDER::Scatterer::get_Qback)
      .def_property_readonly("Qpr",   &CYLINDER::Scatterer::get_Qpr)

      .def_property_readonly("Csca",  &CYLINDER::Scatterer::get_Csca)
      .def_property_readonly("Cext",  &CYLINDER::Scatterer::get_Cext)
      .def_property_readonly("Cabs",  &CYLINDER::Scatterer::get_Cabs)
      .def_property_readonly("Cback", &CYLINDER::Scatterer::get_Cback)
      .def_property_readonly("Cpr",   &CYLINDER::Scatterer::get_Cpr)

      .def_property_readonly("g",     &CYLINDER::Scatterer::get_g)

      .def_readwrite("state", &CYLINDER::Scatterer::state)
      .def_readwrite("area", &CYLINDER::Scatterer::area)
      .def_readwrite("size_parameter", &CYLINDER::Scatterer::size_parameter)
      ;
}







// -
