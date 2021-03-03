#include "Functions.cpp"
#include <iostream>
namespace py = pybind11;

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::map<int, double> dict;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;
#define j complex128(0.0,1.0)









PYBIND11_MODULE(Sphere, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";


}














//-
