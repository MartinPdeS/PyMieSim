
#include <vector>
#include <complex>
#include <math.h>
#include <tuple>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <complex_bessel.h>

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<double> Vec;

#define JJ std::complex<double>(0.0,1.0)

#define PI (double)3.14159265358979323846264338


double
nmFactorial(int n, int m)
{
  double mtemp = boost::math::lgamma(n+m+1);
  double ntemp = boost::math::lgamma(n-m+1);
  return exp(ntemp - mtemp);
}

double Pnm(int n, int m, double x){return boost::math::legendre_p(n, m, x); }
double NPnm(int n, int m, double x){return sqrt((2.*(double)n + 1.)/2. * nmFactorial(n,m)) * Pnm(n,m,x); }
double Pnm_p(int n, int m, double x){ return (sqrt(1-x*x) * Pnm(n,m+1,x) + m*x*Pnm(n,m,x))/(x*x-1); }

double
Pinm(int n, int m, double x)
{
  if (x >= 1-1e-6){x = 1-1e-6;}
  if (x <= -1+1e-6){x = -1+1e-6;}
  return -Pnm(n,m,x) / sqrt(1-x*x);
}


double
Taunm(int n, int m, double x)
{
  if (x >= 1-1e-6){x = 1-1e-6;}
  if (x <= -1+1e-6){x = -1+1e-6;}
  return sqrt(1-x*x) * Pnm_p(n, m, x);
}


//---------------------------------complex-arg-bessel--------------------------------------
template<typename T>
inline complex128 F90jn(int order, T x){ return sp_bessel::sph_besselJ(order, x); }
template<typename T>
inline complex128 F90yn(int order, T x){ return sp_bessel::sph_besselY(order, x); }

template<typename T>
inline complex128 F90Jn(int order, T x){ return sp_bessel::besselJ(order, x); }
template<typename T>
inline complex128 F90Yn(int order, T x){ return sp_bessel::besselY(order, x); }

template<typename T>
inline complex128 F90Jn_p(int order, T x){ return sp_bessel::besselJp(order, x, 1); }
template<typename T>
inline complex128 F90Yn_p(int order, T x){ return sp_bessel::besselYp(order, x, 1); }

template<typename T>
inline complex128 F90h1(int order, T x){ return sp_bessel::sph_hankelH1(order, x); }
template<typename T>
inline complex128 F90h2(int order, T x){ return sp_bessel::sph_hankelH2(order, x); }

template<typename T>
inline complex128 F90H1(int order, T x){ return sp_bessel::hankelH1(order, x); }
template<typename T>
inline complex128 F90H2(int order, T x){ return sp_bessel::hankelH2(order, x); }

template<typename T>
inline complex128 F90H1_p(int order, T x){ return sp_bessel::hankelH1p(order, x); }
template<typename T>
inline complex128 F90H2_p(int order, T x){ return sp_bessel::hankelH2p(order, x); }

//---------------------------------real-arg-bessel--------------------------------------
inline double jn(int order, double x){ return boost::math::sph_bessel(order, x); }
inline double yn(int order, double x){ return boost::math::sph_neumann(order, x); }

inline double yn_p(int order, double x){ return boost::math::sph_neumann_prime(order, x); }
inline double jn_p(int order, double x){ return boost::math::sph_bessel_prime(order, x); }

inline double jn(double order, double x){ return boost::math::sph_bessel(order, x); }
inline double yn(double order, double x){ return boost::math::sph_neumann(order, x); }

inline double yn_p(double order, double x){ return boost::math::sph_neumann_prime(order, x); }
inline double jn_p(double order, double x){ return boost::math::sph_bessel_prime(order, x); }

inline double  Yn(int order, double x){ return boost::math::cyl_neumann(order, x); }
inline double  Jn(int order, double x){ return boost::math::cyl_bessel_j(order, x); }

inline double Yn_p(int order, double x){ return boost::math::cyl_neumann_prime(order, x); }
inline double Jn_p(int order, double x){ return boost::math::cyl_bessel_j_prime(order, x); }

inline double Yn(double order, double x){ return boost::math::cyl_neumann(order, x); }
inline double Jn(double order, double x){ return boost::math::cyl_bessel_j(order, x); }

inline double Yn_p(double order, double x){ return boost::math::cyl_neumann_prime(order, x); }
inline double Jn_p(double order, double x){ return boost::math::cyl_bessel_j_prime(order, x); }

inline complex128 Hn(int order, double x){ return Jn(order,x) + JJ * Yn(order,x); }
inline complex128 Hn_p(int order, double x){ return Jn_p(order,x) + JJ * Yn_p(order,x); }

std::pair<std::vector<double> , std::vector<double>>
Arrange(const double start, const double stop)
{
  std::vector<double> Vec0 ;
  std::vector<double> Vec1 ;
  for (double i = start; i < stop; i++)
  {
    Vec0.push_back(i);
    Vec1.push_back( ( 2 * (i) + 1) / ( (i) * (i + 1) ) ) ;
  }
  return std::make_pair(Vec0, Vec1);
}

inline
complex128 _Psi_p(int type, int n, double x)
{
  if (type == 0){return (complex128) (x * jn_p(n, x) + jn(n, x));}
  if (type == 1){return (complex128) jn_p(n, x);}
  if (type == 2){return (complex128) yn_p(n, x);}
  if (type == 3){return (complex128) jn_p(n, x) + JJ * yn_p(n, x);}
  if (type == 4){return (complex128) jn_p(n, x) - JJ * yn_p(n, x);}

  return 0.;
}

inline
complex128 _Psi(int type, int n, double x)
{
  if (type == 0){return (complex128) x * jn(n, x) ;}
  if (type == 1){return (complex128) jn(n, x) ;}
  if (type == 2){return (complex128) yn(n, x) ;}
  if (type == 3){return (complex128) _Psi(1, n, x) + JJ * _Psi(2, n, x) ;}
  if (type == 4){return (complex128) _Psi(1, n, x) - JJ * _Psi(2, n, x) ;}

  return 0.;
}


inline complex128 Psi(int n, double x){ return (complex128) x * _Psi(1, n, x) ; }
inline complex128 Psi_p(int n, double x){ return (complex128)  x * _Psi_p(1, n, x) + _Psi(1, n, x) ; }

inline complex128 Xi(int n, double x){return x * _Psi(4,n,x); }
inline complex128 Xi_p(int n, double x){return x * _Psi_p(4,n,x) + _Psi(4,n,x); }








template<typename T>
Vec
Linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);

  return linspaced;
}

// -
