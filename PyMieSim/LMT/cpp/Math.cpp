
#include <vector>
#include <complex>
#include <tuple>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<double> Vec;


#define J std::complex<double>(0.0,1.0)

#define PI (double)3.14159265358979323846264338

int GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }


double jn(int order, double x){ return boost::math::sph_bessel(order, x); }
double yn(int order, double x){ return boost::math::sph_neumann(order, x); }
double yn_p(int order, double x){ return boost::math::sph_neumann_prime(order, x); }
double jn_p(int order, double x){ return boost::math::sph_bessel_prime(order, x); }

double jn(double order, double x){ return boost::math::sph_bessel(order, x); }
double yn(double order, double x){ return boost::math::sph_neumann(order, x); }
double yn_p(double order, double x){ return boost::math::sph_neumann_prime(order, x); }
double jn_p(double order, double x){ return boost::math::sph_bessel_prime(order, x); }


double Yn(int order, double x){ return boost::math::cyl_neumann(order, x); }
double Jn(int order, double x){ return boost::math::cyl_bessel_j(order, x); }
double Yn_p(int order, double x){ return boost::math::cyl_neumann_prime(order, x); }
double Jn_p(int order, double x){ return boost::math::cyl_bessel_j_prime(order, x); }

double Yn(double order, double x){ return boost::math::cyl_neumann(order, x); }
double Jn(double order, double x){ return boost::math::cyl_bessel_j(order, x); }
double Yn_p(double order, double x){ return boost::math::cyl_neumann_prime(order, x); }
double Jn_p(double order, double x){ return boost::math::cyl_bessel_j_prime(order, x); }




complex128 Hn(int order, double x){ return Jn(order,x) + complex128(0.0,1.0) * Yn(order,x); }
complex128 Hn_p(int order, double x){ return Jn_p(order,x) + complex128(0.0,1.0) * Yn_p(order,x); }





std::pair<std::vector<double> , std::vector<double>>
Arrange(const double start,
        const double stop)
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



complex128 _Psi_p(int type, int n, double x)
{
  if (type == 0){return (complex128) (x * jn_p(n, x) + jn(n, x));}
  if (type == 1){return (complex128) jn_p(n, x);}
  if (type == 2){return (complex128) yn_p(n, x);}
  if (type == 3){return (complex128) jn_p(n, x) + J * yn_p(n, x);}
  if (type == 4){return (complex128) jn_p(n, x) - J * yn_p(n, x);}

  return 0.;
}

complex128 _Psi(int type, int n, double x)
{
  if (type == 0){return (complex128) x * jn(n, x) ;}
  if (type == 1){return (complex128) jn(n, x) ;}
  if (type == 2){return (complex128) yn(n, x) ;}
  if (type == 3){return (complex128) _Psi(1, n, x) + J * _Psi(2, n, x) ;}
  if (type == 4){return (complex128) _Psi(1, n, x) - J * _Psi(2, n, x) ;}

  return 0.;
}


complex128 Psi(int n, double x){ return (complex128) x * _Psi(1, n, x) ; }

complex128 Psi_p(int n, double x){ return (complex128)  x * _Psi_p(1, n, x) + _Psi(1, n, x) ; }


complex128 Xi(int n, double x){return x * _Psi(4,n,x); }

complex128 Xi_p(int n, double x){return x * _Psi_p(4,n,x) + _Psi(4,n,x); }






// -
