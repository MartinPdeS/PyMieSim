#pragma once

#include <vector>
#include <complex>
#include <math.h>
#include <tuple>
#include "../../../libraries/complex_bessel.cpp"

typedef std::complex<double> complex128;

//---------------------------------AMOS_LIBRARY_WRAPPING--------------------------------------
template<typename T, typename U> inline complex128 compute_jn(U order, T x){ return sp_bessel::sph_besselJ(order, x); }
template<typename T, typename U> inline complex128 compute_yn(U order, T x){ return sp_bessel::sph_besselY(order, x); }

template<typename T, typename U> inline complex128 compute_jn_p(U order, T x)  // https://dlmf.nist.gov/10.51 
{
  return sp_bessel::sph_besselJ(order-1, x) - (order+1.0)/x * sp_bessel::sph_besselJ(order, x);
}

template<typename T, typename U> inline complex128 compute_yn_p(U order, T x)
{
  return sp_bessel::sph_besselY(order-1, x) - (order+1.0)/x * sp_bessel::sph_besselY(order, x);
 }


template<typename T, typename U> inline complex128 compute_h1_p(U order, T x) // https://dlmf.nist.gov/10.51
{
 return sp_bessel::sph_hankelH1(order-1, x) - (order+1.0)/x * sp_bessel::sph_hankelH1(order, x);
}

template<typename T, typename U> inline complex128 compute_h2_p(U order, T x)
{
 return sp_bessel::sph_hankelH2(order-1, x) - (order+1.0)/x * sp_bessel::sph_hankelH2(order, x);
}


template<typename T, typename U> inline complex128 compute_Jn(U order, T x){ return sp_bessel::besselJ(order, x); }
template<typename T, typename U> inline complex128 compute_Yn(U order, T x){ return sp_bessel::besselY(order, x); }

template<typename T> inline complex128 compute_Jn_p(int order, T x){ return sp_bessel::besselJp(order, x, 1); }
template<typename T> inline complex128 compute_Yn_p(int order, T x){ return sp_bessel::besselYp(order, x, 1); }

template<typename T> inline complex128 compute_h1(int order, T x){ return sp_bessel::sph_hankelH1(order, x); }
template<typename T> inline complex128 compute_h2(int order, T x){ return sp_bessel::sph_hankelH2(order, x); }

template<typename T> inline complex128 compute_H1(int order, T x){ return sp_bessel::hankelH1(order, x); }
template<typename T> inline complex128 compute_H2(int order, T x){ return sp_bessel::hankelH2(order, x); }

template<typename T> inline complex128 compute_H1_p(int order, T x){ return sp_bessel::hankelH1p(order, x); }
template<typename T> inline complex128 compute_H2_p(int order, T x){ return sp_bessel::hankelH2p(order, x); }


// -
