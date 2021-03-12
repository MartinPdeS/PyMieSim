#include <vector>
#include <complex>
#include <boost/math/special_functions.hpp>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;









// -
