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

using namespace std;

#define j complex128(0.0,1.0)

int GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }



void
CoefficientAnBn(const double &SizeParam,
                const double &Index,
                const double &nMedium,
                const int    &MaxOrder,
                complex128   *an,
                complex128   *bn)
{
  double alpha = SizeParam,
         beta  = alpha * Index,
         MuSp  = 1.,
         Mu    = 1.,
         M     = Index/nMedium;

  complex128 numerator, denominator, PsiAlpha, PsiBeta, PsiPBeta, PsiPAlpha, XiAlpha, XiPAlpha;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    PsiAlpha  = Psi(order, alpha);
    PsiBeta   = Psi(order, beta);
    PsiPBeta  = Psi_p(order, beta);
    PsiPAlpha = Psi_p(order, alpha);
    XiAlpha   = Xi(order, alpha);
    XiPAlpha  = Xi_p(order, alpha);

    numerator = MuSp * PsiAlpha * PsiPBeta  - Mu * M * PsiPAlpha * PsiBeta;
    denominator = MuSp * XiAlpha * PsiPBeta - Mu * M * XiPAlpha * PsiBeta;
    an[order-1] = numerator/denominator;

    numerator = Mu * M * PsiAlpha * PsiPBeta - MuSp * PsiPAlpha * PsiBeta;
    denominator = Mu * M * XiAlpha * PsiPBeta - MuSp  * XiPAlpha * PsiBeta;
    bn[order-1] = numerator/denominator;
  }
}


iVec
_an(double SizeParam, double Index, double nMedium, int MaxOrder)
{
  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  iVec _an;
  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = MuSp * Psi(order, alpha) * Psi_p(order, beta)  - Mu * M * Psi_p(order, alpha) * Psi(order, beta);
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    _an.push_back(numerator/denominator);
  }
  return _an;
}


iVec
_bn(double SizeParam, double Index, double nMedium, int MaxOrder)
{
  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  iVec _bn;
  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M * Psi(order, alpha) * Psi_p(order, beta) - MuSp * Psi_p(order, alpha) * Psi(order, beta);
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp  * Xi_p(order, alpha) * Psi(order, beta);
    _bn.push_back(numerator/denominator);
  }
  return _bn;
}



iVec
_cn(double SizeParam, double Index, double nMedium, int MaxOrder)
{
  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  iVec _cn;
  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = M * MuSp * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    _cn.push_back(numerator/denominator);
  }
  return _cn;
}


iVec
_dn(double SizeParam, double Index, double nMedium, int MaxOrder)
{

  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  iVec _dn;
  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M*M * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp * Xi_p(order, alpha) * Psi(order, beta);
    _dn.push_back(numerator/denominator);
  }
  return _dn;
}



/*
std::tuple<complex128,complex128>
Expansion(Cndarray   BSC,
          complex128 *_an,
          complex128 *_bn,
          double     Phi,
          double     Theta)
{

  return std::make_tuple(1., 1.);
  py::buffer_info BSCBuffer = BSC.request();


  complex128 *BSCPtr = (complex128 *) BSCBuffer.ptr,
             S1=0.,
             S2=0.,
             TE,
             TM,
             _exp,
             __an,
             __bn;

  double prefactor,
         n_f,
         m_f,
         mu = cos(Phi);

  int Length = BSCBuffer.shape[0],
      Width = BSCBuffer.shape[1],
      m = 0,
      n = 0,
      first = (int)BSCPtr[0].real();

  for (auto i = 0; i < Length; i++)
  {

      n  = (int)BSCPtr[i].real();

      m  = (int)BSCPtr[i + Length].real();

      TE = BSCPtr[i + Length*2];
      TM = BSCPtr[i + Length*3];
      n_f = (double)n; m_f = (double)m;

      prefactor = (2.*n_f+1.)/( n_f* (n_f+1.) );

      _exp   = exp(j*m_f*Theta);
      __an   = _an[n-first+1];
      __bn   = _bn[n-first+1];


      S1 += prefactor * ( ( j * __bn * TE * Taunm[i]) + (m_f * __an * TM * Pinm[i] ) ) * _exp;

      S2 += prefactor * ( ( j * m_f * __bn * TE * Pinm[i]) + (__an * TM * Taunm[i] ) ) * _exp;


  }
  return std::make_tuple(S1, S2);
}


*/





















//-
