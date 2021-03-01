#include "Functions.cpp"
namespace py = pybind11;

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::map<int, double> dict;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;
#define j complex128(0.0,1.0)

Cndarray
S1(double   Index,
   double   Diameter,
   double   Wavelength,
   double   nMedium,
   ndarray  Phi,
   ndarray  Theta,
   double   Polarization,
   double   E0,
   double   R,
   int      Lenght,
   Cndarray  BSC)
{
  py::buffer_info PhiBuffer     = Phi.request();
  py::buffer_info ThetaBuffer   = Theta.request();
  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
         SizeParam = PI * Diameter / Wavelength;

  int PhiLenght   = PhiBuffer.size,
      ThetaLenght = ThetaBuffer.size,
      TotalLenght = PhiLenght*ThetaLenght,
      MaxOrder = GetMaxOrder(SizeParam);

  iVec vec,
      _an = an(SizeParam, Index, nMedium),
      _bn = bn(SizeParam, Index, nMedium);

  Cndarray s1 = ndarray(TotalLenght,0),
           s2 = ndarray(TotalLenght,0);

  auto s1_data = s1.mutable_data();

  int p = 0;
  for(auto i = 0; i < ThetaLenght; i++)
    {
      for(auto l = 0; l < PhiLenght; l++)
        {
          s1_data[p] =  Expansion(MaxOrder, BSC, _an, _bn, PhiPtr[l], ThetaPtr[i]) ;

          p++;
        }
    }

  return s1;
}

 



PYBIND11_MODULE(Sphere, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";



    module.def("S1",
               &S1,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("Polarization"),
               py::arg("E0"),
               py::arg("R"),
               py::arg("Lenght"),
               py::arg("BSC"),
               "Return S1");

}

/*

module.def("an",
           &an,
           py::arg("SizeParam"),
           py::arg("Index"),
           py::arg("nMedium"),
           "Return an");

module.def("bn",
           &bn,
           py::arg("SizeParam"),
           py::arg("Index"),
           py::arg("nMedium"),
           "Return bn");

module.def("cn",
           &cn,
           py::arg("SizeParam"),
           py::arg("Index"),
           py::arg("nMedium"),
           "Return cn");

module.def("dn",
           &dn,
           py::arg("SizeParam"),
           py::arg("Index"),
           py::arg("nMedium"),
           "Return dn");



*/













//-
