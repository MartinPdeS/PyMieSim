#include <boost/math/special_functions/bessel_prime.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#include <cmath>
#include <tuple>




std::tuple<iVec, iVec>
SphereCoefficient(double Radius,
                  std::size_t Order,
                  double Eps,
                  double Mu,
                  double Wavelength)
{

}












PYBIND11_MODULE(GLMT, module) {
    py::module Sphere = module.def_submodule("Sphere", "Sphere scattering object");
    Sphere.def("Coefficient", &SphereCoefficient, "Return coefficient an & bn");
    Sphere.def("ScatteredField", &ScatteredField, "Return scattered field");

}
















//-
