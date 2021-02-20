#include <boost/math/special_functions/bessel_prime.hpp>

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
