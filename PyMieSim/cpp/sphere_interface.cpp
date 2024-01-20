#include <pybind11/pybind11.h>

#include "sphere.cpp"

PYBIND11_MODULE(SphereInterface, module)
{
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";


      pybind11::class_<SPHERE::State>(module, "CppSphereState");


      pybind11::class_<SPHERE::Scatterer>(module, "SPHERE")
      .def(
            pybind11::init<double&, double&, double&, complex128&, double&, CVector&>(),
            pybind11::arg("wavelength"),
            pybind11::arg("amplitude"),
            pybind11::arg("diameter"),
            pybind11::arg("index"),
            pybind11::arg("n_medium"),
            pybind11::arg("jones_vector")
        )

      .def(pybind11::init<>())

      .def(
        "get_s1s2",
        &SPHERE::Scatterer::get_s1s2_py,
        pybind11::arg("phi")
        )

      .def(
            "get_fields",
            &SPHERE::Scatterer::get_unstructured_fields_py,
            pybind11::arg("phi"),
            pybind11::arg("theta"),
            pybind11::arg("r")
        )

      .def(
            "get_full_fields",
            &SPHERE::Scatterer::get_full_structured_fields_py,
            pybind11::arg("sampling"),
            pybind11::arg("r") 
        )


      .def("an", pybind11::overload_cast<>(&SPHERE::Scatterer::get_an_py))
      .def("bn", pybind11::overload_cast<>(&SPHERE::Scatterer::get_bn_py))
      .def("cn", pybind11::overload_cast<>(&SPHERE::Scatterer::get_cn_py))
      .def("dn", pybind11::overload_cast<>(&SPHERE::Scatterer::get_dn_py))

      .def_property_readonly("Qsca",  &SPHERE::Scatterer::get_Qsca)
      .def_property_readonly("Qext",  &SPHERE::Scatterer::get_Qext)
      .def_property_readonly("Qabs",  &SPHERE::Scatterer::get_Qabs)
      .def_property_readonly("Qback", &SPHERE::Scatterer::get_Qback)
      .def_property_readonly("Qpr",   &SPHERE::Scatterer::get_Qpr)


      .def_property_readonly("Csca",  &SPHERE::Scatterer::get_Csca)
      .def_property_readonly("Cext",  &SPHERE::Scatterer::get_Cext)
      .def_property_readonly("Cabs",  &SPHERE::Scatterer::get_Cabs)
      .def_property_readonly("Cback", &SPHERE::Scatterer::get_Cback)
      .def_property_readonly("Cpr",   &SPHERE::Scatterer::get_Cpr)

      .def_property_readonly("g", &SPHERE::Scatterer::get_g)


      .def_readwrite("state", &SPHERE::Scatterer::state)
      .def_readwrite("area", &SPHERE::Scatterer::area)
      .def_readwrite("size_parameter", &SPHERE::Scatterer::size_parameter)
      ;
}







// -
