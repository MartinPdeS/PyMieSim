#include <boost/math/special_functions/bessel_prime.hpp>
#include <eigen3/Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#include "Math.cpp"
#include "Special.cpp"
#include "utils.cpp"
#include <cmath>
#include <tuple>

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<std::vector<double>> Matrix;
typedef std::map<int, double> dict;

int GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }


iVec*
an(double SizeParam, double Index, double MediumIndex)
{
  int MaxOrder = GetMaxOrder(SizeParam);
  double alpha = SizeParam, beta = alpha * Index;
  double MuSp = 1., Mu = 1., M = Index/MediumIndex;
  complex128 numerator, denominator;
  iVec* _an = new iVec(MaxOrder);
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = MuSp * Psi(order, alpha) * Psi_p(order, beta)  - Mu * M * Psi_p(order, alpha) * Psi(order, beta);
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    _an->push_back(numerator/denominator);
  }
  return _an;
}


iVec*
bn(double SizeParam, double Index, double MediumIndex)
{
  int MaxOrder = GetMaxOrder(SizeParam);
  double alpha = SizeParam, beta = alpha * Index;
  double MuSp = 1., Mu = 1., M = Index/MediumIndex;
  complex128 numerator, denominator;
  iVec* _bn = new iVec(MaxOrder);
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M * Psi(order, alpha) * Psi_p(order, beta) - MuSp * Psi_p(order, alpha) * Psi(order, beta);
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp  * Xi_p(order, alpha) * Psi(order, beta);
    _bn->push_back(numerator/denominator);
  }
  return _bn;
}



iVec*
cn(double SizeParam, double Index, double MediumIndex)
{
  int MaxOrder = GetMaxOrder(SizeParam);
  double alpha = SizeParam, beta = alpha * Index;
  double MuSp = 1., Mu = 1., M = Index/MediumIndex;
  complex128 numerator, denominator;
  iVec* _cn = new iVec(MaxOrder);
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = M * MuSp * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    _cn->push_back(numerator/denominator);
  }
  return _cn;
}


iVec*
dn(double SizeParam, double Index, double MediumIndex)
{
  int MaxOrder = GetMaxOrder(SizeParam);
  double alpha = SizeParam, beta = alpha * Index;
  double MuSp = 1., Mu = 1., M = Index/MediumIndex;
  complex128 numerator, denominator;
  iVec* _dn = new iVec(MaxOrder);
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M*M * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp * Xi_p(order, alpha) * Psi(order, beta);
    _dn->push_back(numerator/denominator);
  }
  return _dn;
}


complex128
Expansion(int    MaxOrder,
          dict   BSC_TE,
          dict   BSC_TM,
          iVec*  _an,
          iVec*  _bn,
          double Phi,
          double Theta)
{
  complex128 j = (complex128) 1j;
  double prefactor, TE, TM, order_d;

  complex128 Lterm, Rterm, result=0.;
  for (auto order = 1; order < MaxOrder+1; order++)
  {
      order_d = (double)order;
      prefactor = (2.*order_d+1.)/( order_d* (order_d+1.) );

      for (auto m = -order; m < order+1; m++)
      {
        TE = BSC_TE[m]; TM = BSC_TM[m];
        if (TE == 0. || TM == 0.){continue;}
        Lterm = order_d * (*_an)[order] * TM * Pinm(order, abs(m), cos(Phi - PI/2));
        Rterm = j * (*_bn)[order] * TE * Taunm(order, abs(m), cos(Phi - PI/2));
        result += prefactor*(Rterm+Lterm) * exp(j*order_d*Theta);
      }
  }
  return result;
}



iVec
S1(double   SizeParam,
   double   Index,
   double   MediumIndex,
   Vec      PhiList,
   Vec      ThetaList,
   dict     BSC_TM,
   dict     BSC_TE)
{

  iVec vec;
  complex128 j = (complex128) 1j;


  int MaxOrder = GetMaxOrder(SizeParam);

  iVec* _an = an(SizeParam, Index, MediumIndex);
  iVec* _bn = bn(SizeParam, Index, MediumIndex);

  for(auto const& Theta: ThetaList)
    {
      for(auto const& Phi: PhiList)
        {
          vec.push_back( Expansion(MaxOrder, BSC_TE, BSC_TM, _an, _bn, Phi, Theta) );
        }
    }

  return vec;
}





PYBIND11_MODULE(Sphere, module) {
    //py::module Sphere = module.def_submodule("Sphere", "Sphere scattering object");
    module.def("an", &an, "Return an");
    module.def("bn", &bn, "Return bn");
    module.def("cn", &cn, "Return cn");
    module.def("dn", &dn, "Return dn");
    module.def("S1", &S1, "Return S1");

}















//-
