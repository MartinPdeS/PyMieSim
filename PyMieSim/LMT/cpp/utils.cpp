#include <vector>
#include <complex>
#include <cmath>


typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<double> Vec;

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;
typedef py::buffer_info info;
typedef unsigned int uint;

uint
GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }






double
GetSizeParameter(const double Diameter,
                 const double Wavelength,
                 const double nMedium)
{
  return PI * Diameter / (Wavelength / nMedium);
}
