#include <vector>
#include <complex>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

typedef py::array_t<double> ndarray;
typedef py::array_t<std::complex<double>> Cndarray;
typedef py::buffer_info info;


uint
GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }

double
GetSizeParameter(const double Diameter,
                 const double Wavelength,
                 const double nMedium)
{
  return PI * Diameter / (Wavelength / nMedium);
}



template <class T>
T Sum(const std::vector<T>* vector)
{
   const long unsigned int N = vector->size();
   T sum = 0.;
   for (long unsigned int i = 0; i < N; i++)
   {
     sum += (*vector)[i];
   }
   return sum;
}


std::vector<double>
linespace(const double start,
          const double end,
          const long unsigned int N)
{
    std::vector<double> vector = std::vector<double>(N);

    const double delta = (end-start)/N;

    for (long unsigned int i = 0; i < N; i++)
      {
        vector[i] = delta*i;
      }
      return vector;
}
