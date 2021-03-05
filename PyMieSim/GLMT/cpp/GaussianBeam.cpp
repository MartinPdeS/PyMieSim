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

#define EPS 1e-10

using boost::math::quadrature::trapezoidal;
using boost::math::constants::two_pi;
using boost::math::constants::pi;


struct argument
  {
  int n, m;
  double rhon, angle, s, R0;
  complex128 Q;
  Vec Offset;
} args, args0 ;



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

 complex128 integral = simpson_rule(func,  0.0, two_pi<double>(), 1000);

 return 1./( two_pi<double>()) * integral;
}

complex128
ImTrapz(int m, complex128 beta)
{
 auto func = [=](double angle){return exp( beta * cos(angle)  - j *  (double)m * angle) ;};

 complex128 integral = trapezoidal(func, 0.0, two_pi<double>(),EPS);

 return 1./( two_pi<double>()) * integral;
}



complex128
ImHat(double     m,
      complex128 beta,
      double     xi)
{
  return ImTrapz(m, beta) * exp(-beta - j * m * xi ) ;
}




complex128
I_0(argument arg)
{

  complex128 term0 = pow(-j, arg.n) * arg.rhon*arg.rhon,
             term1 = (2. * (double)arg.n + 1.) * Xi(arg.n, arg.rhon),
             term2 = exp(-j * arg.Offset[2]);

  return term0 / term1 * term2;
}

complex128 I_2(argument arg){ return (2. * arg.Q * arg.s*arg.s * arg.rhon * cos(arg.angle) - 1. ) * sin(arg.angle); }

complex128 I_3(argument arg){ return 4. * arg.Q * arg.s*arg.s * arg.Offset[0] * cos(arg.angle); }

complex128 I_4(argument arg){ return 4. * arg.Q * arg.s*arg.s * arg.Offset[1] * cos(arg.angle); }

complex128
I_1(argument arg)
{
  complex128 term0 = j * arg.Q * arg.s*arg.s * pow(arg.R0 - arg.rhon * sin(arg.angle), 2. ),
             term2 = j * arg.rhon * cos(arg.angle),
             term3 = NPnm(arg.n, abs(arg.m), cos(arg.angle)) * sin(arg.angle) ;

  return arg.Q * exp(term0 + term2 ) * term3;
}


complex128
Anm_integrand(double angle,
              int n,
              int m,
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

  _Q   = Q(r, angle, w0, offset, k);

  args.angle = angle; args.n=n; args.m=m; args.rhon=rhon; args.Offset=Offset; args.s=s; args.Q=_Q; args.R0=R0;

  beta = -2. * j * _Q * s*s * R0 * rhon * sin(angle);

  term0 =  ImHat(m+1, beta, xi);

  term0 += ImHat(m-1, beta, xi);

  term0 *= -I_2(args);

  term0 -= (I_3(args) * ImHat(m, beta, xi));

  term0 *= I_1(args);

  return term0 ;

}



complex128
Bnm_integrand(double angle,
              int n,
              int m,
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

  _Q   = Q(r, angle, w0, offset, k);

  args.angle = angle; args.n=n; args.m=m; args.rhon=rhon; args.Offset=Offset; args.s=s; args.Q=_Q; args.R0=R0;

  beta = -2. * j * _Q * s*s * R0 * rhon * sin(angle);

  term0 =  ImHat(m+1, beta, xi);

  term0 -= ImHat(m-1, beta, xi);

  term0 *= -I_2(args);

  term0 -= (I_4(args) * ImHat(m, beta, xi));

  term0 *= I_1(args);

  return term0 ;

}
 

complex128
Anm(int n,
    int m,
    double k,
    double w0,
    double s,
    Vec Offset,
    Vec offset,
    double R0,
    double xi)
{
  double rhon = (n + 0.5);
  argument args0; args0.n = n; args0.rhon = rhon; args0.Offset = Offset;

  auto func_ = [=](double angle) {return Anm_integrand(angle, n, m, k, w0, s, Offset, offset, R0, xi);};

  return -trapezoidal(func_, 0.0, pi<double>(),EPS) * I_0(args0);

}



complex128
Bnm(int n,
    int m,
    double k,
    double w0,
    double s,
    Vec Offset,
    Vec offset,
    double R0,
    double xi)
{
  double rhon = (n + 0.5);
  argument args0; args0.n = n; args0.rhon = rhon; args0.Offset = Offset;

  auto func_ = [=](double angle) {return Bnm_integrand(angle, n, m, k, w0, s, Offset, offset, R0, xi);};

  return -trapezoidal(func_, 0.0, pi<double>(),EPS) * I_0(args0);

}



std::tuple<ndarray,Cndarray>
PyAnm_integrand(int n,
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

  Cndarray Angle = Cndarray(sampling);
  auto Angle_data = Angle.mutable_data();

  Vec X = Linspace( 0.0, pi<double>(),sampling);

  for (auto i=0; i<sampling;i++)
  {
    _Anm_data[i] = Anm_integrand(X[i], n, m, k, w0, s, Offset, offset, R0, xi);
    Angle_data[i] = X[i];
  }
  return std::make_tuple(Angle, _Anm);
}


std::tuple<ndarray,Cndarray>
PyBnm_integrand(int n,
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

  Cndarray Angle = Cndarray(sampling);
  auto Angle_data = Angle.mutable_data();

  Vec X = Linspace( 0.0, pi<double>(),sampling);

  for (auto i=0; i<sampling;i++)
  {
    _Anm_data[i] = Bnm_integrand(X[i], n, m, k, w0, s, Offset, offset, R0, xi);
    Angle_data[i] = X[i];
  }
  return std::make_tuple(Angle, _Anm);
}





PYBIND11_MODULE(GaussianBeam, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

    module.def("Anm",
               &Anm,
               py::arg("n"),
               py::arg("m"),
               py::arg("k"),
               py::arg("w0"),
               py::arg("s"),
               py::arg("Offset"),
               py::arg("offset"),
               py::arg("R0"),
               py::arg("xi"),
               "Compute Anm");


     module.def("Anm_integrand",
                &PyAnm_integrand,
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
                "Compute Anm integrand as a function of theta");



    module.def("Bnm",
               &Bnm,
               py::arg("n"),
               py::arg("m"),
               py::arg("k"),
               py::arg("w0"),
               py::arg("s"),
               py::arg("Offset"),
               py::arg("offset"),
               py::arg("R0"),
               py::arg("xi"),
               "Compute Bnm");


     module.def("Bnm_integrand",
                &PyBnm_integrand,
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
                "Compute Bnm integrand as a function of theta");
}














//-
