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










iVec
an(double SizeParam, double Index, double nMedium)
{
  int MaxOrder = GetMaxOrder(SizeParam);
  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium; 

  complex128 numerator, denominator;
  iVec _an;
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = MuSp * Psi(order, alpha) * Psi_p(order, beta)  - Mu * M * Psi_p(order, alpha) * Psi(order, beta);
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    _an.push_back(numerator/denominator);

  }
  return _an;
}


iVec
bn(double SizeParam, double Index, double nMedium)
{
  int MaxOrder = GetMaxOrder(SizeParam);

  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  iVec _bn;
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M * Psi(order, alpha) * Psi_p(order, beta) - MuSp * Psi_p(order, alpha) * Psi(order, beta);
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp  * Xi_p(order, alpha) * Psi(order, beta);
    _bn.push_back(numerator/denominator);
  }
  return _bn;
}



iVec
cn(double SizeParam, double Index, double nMedium)
{
  int MaxOrder = GetMaxOrder(SizeParam);

  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  iVec _cn;
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = M * MuSp * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    _cn.push_back(numerator/denominator);
  }
  return _cn;
}


iVec
dn(double SizeParam, double Index, double nMedium)
{
  int MaxOrder = GetMaxOrder(SizeParam);

  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  iVec _dn;
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M*M * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp * Xi_p(order, alpha) * Psi(order, beta);
    _dn.push_back(numerator/denominator);
  }
  return _dn;
}




std::pair<Cndarray, Cndarray>
S1(double   Index,
   double   Diameter,
   double   Wavelength,
   double   nMedium,
   ndarray  Phi,
   ndarray  Theta,
   double   Polarization,
   double   E0,
   double   R,
   Cndarray  BSC,
   int       MaxOrder)
{
  py::buffer_info PhiBuffer     = Phi.request();
  py::buffer_info ThetaBuffer   = Theta.request();
  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
         SizeParam = PI * Diameter / Wavelength;

  int PhiLenght   = PhiBuffer.size,
      ThetaLenght = ThetaBuffer.size,
      TotalLenght = PhiLenght*ThetaLenght;

  iVec vec,
      __an = _an(SizeParam, Index, nMedium, MaxOrder),
      __bn = _bn(SizeParam, Index, nMedium, MaxOrder);

  Cndarray s1 = Cndarray(TotalLenght),
           s2 = Cndarray(TotalLenght);

  complex128 temp0,
             temp1;

  auto s1_data = s1.mutable_data(),
       s2_data = s2.mutable_data();

  int p = 0;
  for(auto i = 0; i < ThetaLenght; i++)
    {
      for(auto l = 0; l < PhiLenght; l++)
        {

          std::tie(temp0, temp1) = Expansion(MaxOrder, BSC, __an, __bn, PhiPtr[l], ThetaPtr[i]) ;

          s1_data[p] = temp0;
          s2_data[p] = temp1;
          p++;
        }
    }

  return std::make_pair(s1, s2);
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
               py::arg("BSC"),
               py::arg("MaxOrder"),
               "Return S1");


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
