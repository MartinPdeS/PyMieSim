#include <boost/math/special_functions/bessel_prime.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#include "Math.cpp"
#include "Special.cpp"
#include "utils.cpp"
#include <cmath>

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::map<int, double> dict;
typedef std::list<std::list<int>> List2D;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;

#define j complex128(0.0,1.0)

int GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }


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
  iVec _an = iVec(MaxOrder);
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
  iVec _bn = iVec(MaxOrder);
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
  iVec _cn = iVec(MaxOrder);
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
  iVec _dn = iVec(MaxOrder);
  for (long unsigned int order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M*M * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp * Xi_p(order, alpha) * Psi(order, beta);
    _dn.push_back(numerator/denominator);
  }
  return _dn;
}


complex128
Expansion(int      MaxOrder,
          Cndarray BSC,
          iVec    _an,
          iVec   _bn,
          double   Phi,
          double   Theta)
{

  py::buffer_info BSCBuffer = BSC.request();


  complex128 *BSCPtr = (complex128 *) BSCBuffer.ptr,
             Lterm,
             Rterm,
             result=0.,
             TE,
             TM;

  double prefactor,
         order_f;

  int Lenght = BSCBuffer.shape[0],
      Width = BSCBuffer.shape[1],
      m = 0,
      n = 0;

  for (auto i = 0; i < Lenght; i++)
  {

      n  = (int)BSCPtr[i*Width + 0].real();
      m  = (int)BSCPtr[i*Width + 1].real();
      TE = BSCPtr[i*Width + 2];
      TM = BSCPtr[i*Width + 3];
      order_f = (double)n;

      prefactor = (2.*order_f+1.)/( order_f* (order_f+1.) );

      //if (TE == 0. || TM == 0.){continue;}
      Lterm = order_f * _an[n] * TM * Pinm(n, abs(m), cos(Phi - PI/2));
      Rterm = j * _bn[n] * TE * Taunm(n, abs(m), cos(Phi - PI/2));
      result += prefactor*(Rterm+Lterm) * exp(j*order_f*Theta);

  }
  return result;
}






/*

for (auto order = 1; order < MaxOrder+1; order++)
{
    order_d = (double)order;
    prefactor = (2.*order_d+1.)/( order_d* (order_d+1.) );

    for (auto m = -order; m < order+1; m++)
    {
      TE = BSCPtr[order*Y  + 2]; TM = BSCPtr[m];
      if (TE == 0. || TM == 0.){continue;}
      Lterm = order_d * _an[order] * TM * Pinm(order, abs(m), cos(Phi - PI/2));
      Rterm = j * _bn[order] * TE * Taunm(order, abs(m), cos(Phi - PI/2));
      result += prefactor*(Rterm+Lterm) * exp(j*order_d*Theta);
    }
}

*/






















//-
