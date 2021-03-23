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

/*
template<class Type, class IndexType>
py::array_t<Type, py::array::c_style> Multi_array<Type,1,IndexType>::move_py(IndexType n_i)
{
	py::capsule free_when_done( ptr, free_func );
	py::array_t<Type, py::array::c_style> out
	(
		{n_i},      // shape
		{stride_i}, // C-style contiguous strides for double
		ptr,     	// the data pointer
		free_when_done // Free method
	);

	alloc_func = NULL ;
	free_func = NULL ;
	ptr = NULL ;
	n_i = 0 ;
	stride_i = 0 ;

	return out ;
};
*/


void
LoopUnstructured(int         PhiLength,
                 int         ThetaLength,
                 complex128  propagator,
                 double     *PhiPtr,
                 double     *ThetaPtr,
                 complex128 *EPhiPtr,
                 complex128 *EThetaPtr,
                 complex128 *s1s2,
                 bool        Polarized,
                 double      Polarization)

{

  double temp0;

  for (auto p=0; p < PhiLength; p++ )
  {
        temp0 = ThetaPtr[p] ;
        EPhiPtr[p] = J* propagator * s1s2[p] * (complex128) abs(cos(temp0 + Polarization));
        EThetaPtr[p] =- propagator * s1s2[p + PhiLength] * (complex128) abs(sin(temp0 + Polarization));
  }
}



void
LoopStructured(int         PhiLength,
               int         ThetaLength,
               complex128  propagator,
               double     *PhiPtr,
               double     *ThetaPtr,
               complex128 *EPhiPtr,
               complex128 *EThetaPtr,
               complex128 *s1s2,
               bool        Polarized,
               double      Polarization)

{
  int w = 0;
  double temp0;

  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          temp0 = ThetaPtr[t] ;
          EPhiPtr[w] = J* propagator * s1s2[p] * (complex128) abs(cos(temp0 + Polarization));
          EThetaPtr[w] =- propagator * s1s2[p + PhiLength] * (complex128) abs(sin(temp0 + Polarization));
          w++;
       }
  }
}




double
GetSizeParameter(const double Diameter,
                 const double Wavelength,
                 const double nMedium)
{
  return PI * Diameter / (Wavelength / nMedium);
}
