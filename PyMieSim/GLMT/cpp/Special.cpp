#include <iostream>
#include <complex>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/special_functions/gamma.hpp>

typedef std::complex<double> complex128;


#define j complex128(0.0,1.0)
#define PI (double) 3.14159265358979323846264338


using namespace std;

std::tuple<Vec, Vec, Vec>
Cart2sph(Vec X, Vec Y, Vec Z){
  Vec R, PHI, THETA;
  for (long unsigned int i = 0; i < X.size(); i++){
    R.push_back( sqrt( pow( X[i],2) + pow(Y[i],2) + pow(Z[i],2) ) );
    PHI.push_back( atan( Y[i]/X[i] ) );
    THETA.push_back( atan( sqrt( pow(Y[i],2) / pow(X[i],2) ) ) );
  }
  return std::make_tuple(R, PHI, THETA);
}

int GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }


std::tuple<Vec, Vec, Vec>
Fibonacci_sphere(int samples,
                 double maxAngle){


    Vec X, Y , Z;

    Vec R, THETA, PHI;
    double phi = PI * (3. - sqrt(5.));


    double SolidAngle = abs( 2*PI * (cos(maxAngle) - 1));

    double ratio = 4*PI / SolidAngle;

    long unsigned int TrueSampling = (long unsigned int)( samples * ratio );

    for (long unsigned int i = 0; i < TrueSampling; i++){

        double y = 1 - (i / (double)(TrueSampling - 1)) * 2 ;
        double radius = sqrt(1 - y * y);

        double theta = phi * i;

        double x = cos(theta) * radius;
        //double z = sin(theta) * radius;

        X.push_back(x);
        Y.push_back(1.0);
        Z.push_back(1.0);

        if (i >= (long unsigned int)samples - 1) break;
      }

    std::tie(R, PHI, THETA) = Cart2sph(X,Y,Z);

    X.clear();
    Y.clear();
    Z.clear();

    X.shrink_to_fit();
    Y.shrink_to_fit();
    Z.shrink_to_fit();

    return std::make_tuple(R, PHI, THETA);


}


double
nmFactorial(int n, int m)
{
  double mtemp = boost::math::lgamma(n+m+1);

  double ntemp = boost::math::lgamma(n-m+1);

  return exp(ntemp - mtemp);
}


double Yn(int order, double x){ return boost::math::cyl_neumann(order, x); }

double Jn(int order, double x){ return boost::math::cyl_bessel_j(order, x); }

double Yn_p(int order, double x){ return boost::math::cyl_neumann_prime(order, x); }

double Jn_p(int order, double x){ return boost::math::cyl_bessel_j_prime(order, x); }

complex128 Hn(int order, double x){ return Jn(order,x) + j * Yn(order,x); }

complex128 Hn_p(int order, double x){ return Jn_p(order,x) + j * Yn_p(order,x); }

double jn(int order, double x){ return boost::math::sph_bessel(order, x); }

double jn_p(int order, double x){ return boost::math::sph_bessel_prime(order, x); }

double yn(int order, double x){ return boost::math::sph_neumann(order, x); }

double yn_p(int order, double x){ return boost::math::sph_neumann_prime(order, x); }

double Pnm(int n, int m, double x){return boost::math::legendre_p(n, m, x); }

double NPnm(int n, int m, double x){return sqrt((2.*(double)n + 1.)/2. * nmFactorial(n,m)) * Pnm(n,m,x); }

complex128 _Psi_p(int type, int n, double x)
{
  if (type == 0){return (complex128) (x * jn_p(n, x) + jn(n, x));}
  if (type == 1){return (complex128) jn_p(n, x);}
  if (type == 2){return (complex128) yn_p(n, x);}
  if (type == 3){return (complex128) jn_p(n, x) + j * yn_p(n, x);}
  if (type == 4){return (complex128) jn_p(n, x) - j * yn_p(n, x);}

  return 0.;
}

complex128 _Psi(int type, int n, double x)
{
  if (type == 0){return (complex128) x * jn(n, x) ;}
  if (type == 1){return (complex128) jn(n, x) ;}
  if (type == 2){return (complex128) yn(n, x) ;}
  if (type == 3){return (complex128) _Psi(1, n, x) + j * _Psi(2, n, x) ;}
  if (type == 4){return (complex128) _Psi(1, n, x) - j * _Psi(2, n, x) ;}

  return 0.;
}


complex128 Psi(int n, double x){ return (complex128) x * _Psi(1, n, x) ; }

complex128 Psi_p(int n, double x){ return (complex128)  x * _Psi_p(1, n, x) + _Psi(1, n, x) ; }


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


complex128 Xi(int n, double x){return x * _Psi(4,n,x); }

complex128 Xi_p(int n, double x){return x * _Psi_p(4,n,x) + _Psi(4,n,x); }








inline void
SymetricPinmTaunm(const int    Length,
                  const double phi,
                  double      *pin,
                  double      *taun)

{

  double mu = cos(phi);

  //if (mu >= 1-1e-6){mu = 1-1e-6;}
  //if (mu <= -1+1e-6){mu = -1+1e-6;}

  double su = sin(phi);


  pin[0] = 0.;

  taun[0] = 0.;

  pin[1] = 1.;

  pin[2] = 3. * mu;

  taun[1] = mu;

  taun[2] = 3. * (mu*mu -su*su );//3. *cos(2*acos(mu));

  double n = 0.;
  for (long unsigned i = 3; i < Length; i++)
      {
       n = (double)i;
       pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

       taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
     }

}









// -
