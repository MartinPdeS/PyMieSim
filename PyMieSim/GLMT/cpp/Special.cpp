
#include <vector>
#include <complex>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<double> Vec;
#define PI 3.14159265


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
Legendre(std::size_t n, double x)
{
  return boost::math::legendre_p(n, 1, x);
}

double
Legendre_p(std::size_t n, double x)
{
  return -boost::math::legendre_p(n, 1, x)/(std::sqrt(1.0 - x*x));
}

double
Legendre_pp(std::size_t n, double x)
{
  return boost::math::legendre_p(n, 2, x)/(1.0 - x*x);
}

double
Pin(std::size_t n, double theta)
{
  return -1 * Legendre_p(n, cos(theta));
}


double
Taun(std::size_t n, double theta)
{
  return - cos(theta) * Legendre_p(n, cos(theta)) + sin(theta)*sin(theta) * Legendre_pp(n, cos(theta));
}



template <class T>
std::vector<T>
Riccati1(std::size_t n, T x)
{
  std::vector<T> temp;
  T b = boost::math::sph_bessel(n, x);
  T b_p = boost::math::sph_bessel_prime(n, x);
  temp.push_back(x * b);
  temp.push_back(b + x * b_p);

  return temp;
}


template <class T>
std::vector<T>
Riccati2(std::size_t n, T x)
{
  std::vector<T> temp;
  T b = boost::math::sph_neumann(n, x);
  T b_p = boost::math::sph_neumann_prime(n, x);
  temp.push_back(x * b);
  temp.push_back(b + x * b_p);

  return temp;
}


iVec
Riccati3(std::size_t n, double x)
{
  complex128 j (0., 1.0);
  iVec ret;
  Vec temp0, temp1;
  temp0 =  Riccati1(n,x);
  temp1 = Riccati2(n,x);


  ret.reserve(temp0.size());

  std::transform( begin(temp0),
                  end(temp0),
                  begin(temp1),
                  std::back_inserter(ret),
                  [](double r, double i) { return std::complex<double>(r, i); }
                 );

  return ret;

}


template <class T>
std::vector<T>
VecHankel(std::size_t n, T x)
{
  complex128 j (0., 1.0);
  iVec ret;
  Vec temp0, temp1;
  temp0 = boost::math::sph_bessel(n, x);
  temp1 = boost::math::sph_neumann(n, x);


  ret.reserve(temp0.size());

  std::transform( begin(temp0),
                  end(temp0),
                  begin(temp1),
                  std::back_inserter(ret),
                  [](double r, double i) { return std::complex<double>(r, i); }
                 );

  return ret;
}


template <class T>
std::vector<T>
VecHankel_p(std::size_t n, T x)
{
  complex128 j (0., 1.0);
  iVec ret;
  Vec temp0, temp1;
  temp0 = boost::math::sph_bessel_prime(n, x);
  temp1 = boost::math::sph_neumann_prime(n, x);


  ret.reserve(temp0.size());

  std::transform( begin(temp0),
                  end(temp0),
                  begin(temp1),
                  std::back_inserter(ret),
                  [](double r, double i) { return std::complex<double>(r, i); }
                 );

  return ret;
}



complex128
Hankel(std::size_t n, double x)
{
  complex128 j (0., 1.0);
  complex128 ret;
  complex128 temp0, temp1;
  temp0 = boost::math::sph_bessel(n, x);
  temp1 = boost::math::sph_neumann(n, x);

  return temp0 + j * temp1;

}


complex128
Hankel_p(std::size_t n, double x)
{
  complex128 j (0., 1.0);
  complex128 ret;
  complex128 temp0, temp1;
  temp0 = boost::math::sph_bessel_prime(n, x);
  temp1 = boost::math::sph_neumann_prime(n, x);

  return temp0 + j * temp1;

}




double jn(int order, double x){ return boost::math::sph_bessel(order, x); }

double jn_p(int order, double x){ return boost::math::sph_bessel_prime(order, x); }

double yn(int order, double x){ return boost::math::sph_neumann(order, x); }

double yn_p(int order, double x){ return boost::math::sph_neumann_prime(order, x); }

double Pnm(int n, int m, double x){return boost::math::legendre_p(n, m, x); }

complex128 _Psi(int type, int n, double x)
{
  if (type == 0){return (complex128) (x * jn_p(n, x) + jn(n, x));}
  if (type == 1){return (complex128) jn_p(n, x);}
  if (type == 2){return (complex128) yn_p(n, x);}
  if (type == 3){return ( jn_p(n, x) + 1.j * yn_p(n,x) ) ;}
  if (type == 4){return ( jn_p(n, x) - 1.j * yn_p(n,x) ) ;}
  return 0.;
}

complex128 _Psi_p(int type, int n, double x)
{
  if (type == 0){return (complex128) x * jn(n, x) ;}
  if (type == 1){return (complex128) jn(n, x) ;}
  if (type == 2){return (complex128) yn(n, x) ;}
  if (type == 4){return (complex128) ( (x * jn_p(n, x) + jn(n, x)) + 1.j * yn_p(n, x) ) ; }
  if (type == 4){return (complex128) ( (x * jn_p(n, x) + jn(n, x)) - 1.j * jn_p(n, x) ) ; }
  return 0.;
}


complex128 Psi(int n, double x){ return (complex128) x * _Psi(1, n, x) ; }

complex128 Psi_p(int n, double x){ return (complex128)  x * _Psi_p(1, n, x) + _Psi(1, n, x) ; }


double Pnm_p(int n, int m, double x){ return (sqrt(1-x*x) * Pnm(n,m+1,x) + m*x*Pnm(n,m,x))/(x*x-1); }


double Pinm(int n, int m, double x)
{
  if (x >= 1-1e-6){x = 1-1e-6;}
  if (x <= -1+1e-6){x = -1+1e-6;}
  return -Pnm(n,m,x) / sqrt(1-x*x);
}


double Taunm(int n, int m, double x)
{
  if (x >= 1-1e-6){x = 1-1e-6;}
  if (x <= -1+1e-6){x = -1+1e-6;}
  return sqrt(1-x*x) * Pnm_p(n, m, x);
}


complex128 Xi(int n, double x){return x * _Psi(4,n,x); }

complex128 Xi_p(int n, double x){return x * _Psi_p(4,n,x) + _Psi(4,n,x); }
































// -
