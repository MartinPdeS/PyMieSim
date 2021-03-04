#include "Functions.cpp"
#include <iostream>
#include <stdio.h>
#include <tgmath.h>
#include <math.h>
#include <complex.h>
namespace py = pybind11;

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::map<int, double> dict;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;
#define j complex128(0.0,1.0)

using boost::math::quadrature::trapezoidal;
using boost::math::constants::two_pi;
using boost::math::constants::pi;


template <typename func_type>
complex128 simpson_rule(func_type f,
                        double a,
                        double b,
                        int n // Number of intervals
                        )
{
    double h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    complex128 sum_odds = 0.0;
    for (int i = 1; i < n; i += 2)
    {
        sum_odds += f(a + i * h);
    }
    complex128 sum_evens = 0.0;
    for (int i = 2; i < n; i += 2)
    {
        sum_evens += f(a + i * h);
    }

    return (f(a) + f(b) + 2. * sum_evens + 4. * sum_odds) * h / 3.;
}



iVec
Q_(double r,
  Vec    theta,
  double w0,
  Vec    offset,
  double k)
{
   iVec result;
   for (auto const& angle: theta)
   {
     result.push_back(1./( 2. * ( r * cos(angle) - offset[2])/(k * w0*w0 ) - j ) ) ;
   }
   return result;
}



complex128
Q(double r,
  double theta,
  double w0,
  Vec    offset,
  double k)
{
  return 1./( 2. * ( r * cos(theta) - offset[2])/(k * w0*w0 ) - j )  ;
}


complex128
ImSimpson(int m, complex128 beta)
{

 auto func = [=](double angle){return exp( beta * cos(angle)  - j *  (double)m * angle) ;};

 complex128 integral = simpson_rule(func,  0.0, two_pi<double>(), 500);

 return 1./( two_pi<double>()) * integral;
}

complex128
ImTrapz(int m, complex128 beta)
{

 auto func = [=](double angle){return exp( beta * cos(angle)  - j *  (double)m * angle) ;};

 complex128 integral = trapezoidal(func, 0.0, two_pi<double>(),1e-5);

 return 1./( two_pi<double>()) * integral;
}


complex128
Im_Sum(int m, complex128 beta)
{

 Vec X = Linspace( 0.0, two_pi<double>(),1000);

 iVec Y;
 double dx = abs(X[1]-X[0]);
 complex128 sum=0.;

 for (auto const& x: X){ sum += exp(-beta*cos(x) - j *  (double)m * x); }

 return sum *dx * 1./( two_pi<double>());

}


complex128
ImHat(double     m,
      complex128 beta,
      double     xi)
{
  return Im_Sum(m, beta) * exp(-beta - j * m * xi ) ;
}




complex128
I_0(int n, int m, double rhon, complex128 Q, Vec Offset)
{
  complex128 term0 = pow(-j, n) * rhon*rhon,
             term1 = (2. * (double)n + 1.) * Xi(n, rhon),
             term2 = exp(-j * Offset[2]);

  return term0 / term1 * term2;
}



complex128
I_2(int n, int m, double rhon, complex128 Q, double theta, Vec Offset, double s, double R0)
{
  return (2. * Q * s*s * rhon * cos(theta) - 1. ) * sin(theta);
}


complex128
I_3(int n, int m, double rhon, complex128 Q, double theta, Vec Offset, double s, double R0)
{
  return 4. * Q * s*s * Offset[0] * cos(theta);
}

complex128
I_4(int n, int m, double rhon, complex128 Q, double theta, Vec Offset, double s, double R0)
{
  return 4. * Q * s*s * Offset[1] * cos(theta);
}

complex128
Integrate(iVec y, Vec Theta)
{
  complex128 sum = 0.;
  double dx  = abs( Theta[1] - Theta[0] );

  for (auto const& val: y)
  {
    sum += val;
  }
  return sum * dx;
}



complex128
I_1(int n, int m, double rhon, complex128 Q, double theta, Vec Offset, double s, double R0)
{
  complex128 term0 = j * Q * s*s,
             term1 = pow(R0 - rhon * sin(theta), 2. ),
             term2 = j * rhon * cos(theta),
             term3 = NPnm(n, abs(m), cos(theta)) * sin(theta) ;

  return Q * exp(term0 * term1 + term2 ) * term3;
}


complex128
Anm_integrand(double angle,
              int n,
              int m,
              int sampling,
              double k,
              double w0,
              double s,
              Vec Offset,
              Vec offset,
              double R0,
              double xi)
{
  if (abs(m) >n){return 0;}

  double rhon = (n + 0.5), r = rhon/k;

  complex128 _Q, beta, term0;
  //beta =2.;angle=2.3;
  _Q   = Q(r, angle, w0, offset, k);

  beta = -2. * j * _Q * s*s * R0 * rhon * sin(angle);

  term0 =  ImHat(m+1, beta, xi);
  term0 += ImHat(m-1, beta, xi);

  term0 *= I_2(n, m, rhon, _Q, angle, Offset, s, R0);

  term0 -= (I_3(n, m, rhon, _Q, angle, Offset, s, R0) * ImHat(m, beta, xi));

  term0 *= I_1(n, m, rhon, _Q, angle, Offset, s, R0);
  //std::cout << term0 << std::endl;
  //std::cout << "term0 is " <<  term0 << std::endl;


  return -term0;

}


complex128
Anm(int n,
    int m,
    int sampling,
    double k,
    double w0,
    double s,
    Vec Offset,
    Vec offset,
    double R0,
    double xi)
{
  //auto func_ = [=](double angle) {return Anm_integrand(angle, n, m, sampling, k, w0, s, Offset, offset, R0, xi);};

  //return trapezoidal(func_, 0.0, pi<double>(),1e-5);
  Vec X = Linspace( 0.0, pi<double>(),400);

  iVec Y;
  double dx = abs(X[1]-X[0]);
  complex128 sum=0.;

  for (auto const& x: X)
  {
    sum += Anm_integrand(x, n, m, sampling, k, w0, s, Offset, offset, R0, xi);
  }
  return sum * dx ;
}



Cndarray
VAnm(int n,
     int m,
     int sampling,
     double k,
     double w0,
     double s,
     Vec Offset,
     Vec offset,
     double R0,
     double xi)
{

  Cndarray _Anm = Cndarray(sampling);
  auto _Anm_data = _Anm.mutable_data();

  Vec X = Linspace( 0.0, pi<double>(),sampling);

  for (auto i=0; i<sampling;i++)
  {
    _Anm_data[i] = Anm_integrand(X[i], n, m, sampling, k, w0, s, Offset, offset, R0, xi);
  }
  return _Anm;
}



PYBIND11_MODULE(GaussianBeam, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

    module.def("Anm",
               &Anm,
               py::arg("n"),
               py::arg("m"),
               py::arg("sampling"),
               py::arg("k"),
               py::arg("w0"),
               py::arg("s"),
               py::arg("Offset"),
               py::arg("offset"),
               py::arg("R0"),
               py::arg("xi"),
               "Compute Anm");


     module.def("VAnm",
                &VAnm,
                py::arg("n"),
                py::arg("m"),
                py::arg("sampling"),
                py::arg("k"),
                py::arg("w0"),
                py::arg("s"),
                py::arg("Offset"),
                py::arg("offset"),
                py::arg("R0"),
                py::arg("xi"),
                "Compute Anm");
}














//-
