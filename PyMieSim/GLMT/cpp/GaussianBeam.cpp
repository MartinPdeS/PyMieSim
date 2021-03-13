#include "Math.cpp"
#include "Special.cpp"
#include <iostream>

#include <stdio.h>
#include <tgmath.h>
#include <math.h>
#include <complex.h>
#include <boost/math/quadrature/trapezoidal.hpp>


#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray; 
typedef py::array_t<complex128> Cndarray;
#define J complex128(0.0,1.0)

#define EPS 1e-10
#define Correction 1.//2.*exp(J*1.0472351410849077)

using boost::math::quadrature::trapezoidal;
using boost::math::constants::two_pi;
using boost::math::constants::pi;


struct argument
  {
  int n, m;
  double r, rhon, angle, s, R0, w0, xi, k;
  complex128 Q, beta;
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
  Vec    Offset,
  double k)
{
   iVec result;
   for (auto const& angle: theta)
   {
     result.push_back(1./( 2. * ( r * cos(angle) - Offset[2]/k )/(k * w0*w0 ) - J ) ) ;
   }
   return result;
}



complex128
Q(double r,
  double theta,
  double w0,
  Vec    Offset,
  double k)
{
  return 1./( 2. * ( r * cos(theta) - Offset[2]/k)/(k * w0*w0 ) - J )  ;
}


complex128
ImSimpson(int m, complex128 beta)
{
 auto func = [=](double angle){return exp( beta * cos(angle)  - J *  (double)m * angle) ;};

 complex128 integral = simpson_rule(func,  0.0, two_pi<double>(), 1000);

 return 1./( two_pi<double>()) * integral;
}

complex128
ImTrapz(int m, complex128 beta)
{
 auto func = [=](double angle)->complex128{return exp( beta * cos(angle)  - J *  (double)m * angle) ;};

 complex128 integral = trapezoidal(func, 0.0, two_pi<double>(), EPS);

 return 1./( two_pi<double>()) * integral;
}



complex128
ImHat(double     m,
      complex128 beta,
      double     xi)
{
  return ImTrapz(m, beta) * exp(-beta - J * m * xi ) ;
}




complex128
I_0(argument arg)
{

  complex128 term0 = pow(-J, arg.n) * arg.rhon*arg.rhon,
             term1 = (2. * (double)arg.n + 1.) * _Psi(0,arg.n, arg.rhon),
             term2 = exp(-J * arg.Offset[2]);

  return term0 / term1 * term2;
}

complex128 I_2(argument arg){ return (2. * arg.Q * arg.s*arg.s * arg.rhon * cos(arg.angle) - 1. ) * sin(arg.angle); }

complex128 I_3(argument arg){ return 4. * arg.Q * arg.s*arg.s * arg.Offset[0] * cos(arg.angle); }

complex128 I_4(argument arg){ return 4. * arg.Q * arg.s*arg.s * arg.Offset[1] * cos(arg.angle); }

complex128
I_1(argument arg)
{
  complex128 term0 = J * arg.Q * arg.s*arg.s * pow(arg.R0 - arg.rhon * sin(arg.angle), 2. ),
             term2 = J * arg.rhon * cos(arg.angle),
             term3 = NPnm(arg.n, abs(arg.m), cos(arg.angle)) * sin(arg.angle) ;

  return arg.Q * exp(term0 + term2 ) * term3;
}






complex128
Bnm_integrand(double angle, argument args0)
{
  if (abs(args0.m) >args0.n){return 0;}

  args0.Q = Q(args0.r, angle, args0.w0, args0.Offset, args0.k);

  args0.angle = angle;

  complex128 beta = -2. * J * args0.Q * args0.s*args0.s * args0.R0 * args0.rhon * sin(angle);

  complex128 term0 =  ImHat(args0.m+1, beta, args0.xi) - ImHat(args0.m-1, beta, args0.xi);

  term0 *= I_2(args0);

  term0 -= (I_4(args0) * ImHat(args0.m, beta, args0.xi));

  term0 *= I_1(args0);

  return term0;

}




complex128
Anm_integrand(double angle, argument args0)
{
  if (abs(args0.m) >args0.n){return 0;}

  args0.Q = Q(args0.r, angle, args0.w0, args0.Offset, args0.k);

  args0.angle = angle;

  complex128 beta = -2. * J * args0.Q * args0.s*args0.s * args0.R0 * args0.rhon * sin(angle);

  complex128 term0 =  ImHat(args0.m+1, beta, args0.xi) + ImHat(args0.m-1, beta, args0.xi);

  term0 *= I_2(args0);

  term0 -= (I_3(args0) * ImHat(args0.m, beta, args0.xi));

  term0 *= I_1(args0);

  return term0;

}



complex128
Anm(int    n,
    int    m,
    double k,
    double w0,
    Vec    Offset,
    double Tolerance)
{
  argument args0;
           args0.n      = n;
           args0.m      = m;
           args0.rhon   = (n + 0.5);
           args0.Offset = Offset;
           args0.k      = k;
           args0.R0     = sqrt(Offset[0]*Offset[0] + Offset[1]*Offset[1]);
           args0.xi     = acos(Offset[0]/args0.R0);
           args0.s      = 1./(k*w0);
           args0.w0     = w0;
           args0.r      = args.rhon/k;

  auto func_ = [=](double angle)->complex128 {return Anm_integrand(angle, args0);};

  return -trapezoidal(func_, 0.0, pi<double>(),Tolerance) * I_0(args0) * Correction;

}



complex128
Bnm(int n,
    int m,
    double k,
    double w0,
    Vec Offset,
    double Tolerance)
{
  argument args0;
           args0.n      = n;
           args0.m      = m;
           args0.rhon   = (n + 0.5);
           args0.Offset = Offset;
           args0.k      = k;
           args0.R0     = sqrt(Offset[0]*Offset[0] + Offset[1]*Offset[1]);
           args0.xi     = acos(Offset[0]/args0.R0);
           args0.s      = 1./(k*w0);
           args0.w0     = w0;
           args0.r      = args.rhon/k;

  auto func_ = [=](double angle)->complex128 {return Bnm_integrand(angle, args0);};

  return -trapezoidal(func_, 0.0, pi<double>(),Tolerance) * I_0(args0) * Correction;

}



std::tuple<ndarray,Cndarray>
PyAnm_integrand(int    n,
                int    m,
                int    sampling,
                double k,
                double w0,
                Vec    Offset,
                double Tolerance)
{

 argument args0;
          args0.n      = n;
          args0.m      = m;
          args0.rhon   = (n + 0.5);
          args0.Offset = Offset;
          args0.k      = k;
          args0.R0     = sqrt(Offset[0]*Offset[0] + Offset[1]*Offset[1]);
          args0.xi     = acos(Offset[0]/args0.R0);
          args0.s      = 1./(k*w0);
          args0.w0     = w0;
          args0.r      = args.rhon/k;

  Cndarray _Anm = Cndarray(sampling);
  auto _Anm_data = _Anm.mutable_data();

  Cndarray Angle = Cndarray(sampling);
  auto Angle_data = Angle.mutable_data();

  Vec X = Linspace( 0.0, pi<double>(),sampling);

  for (auto i=0; i<sampling;i++)
  {
    _Anm_data[i] = Anm_integrand(X[i], args);
    Angle_data[i] = X[i];
  }
  return std::make_tuple(Angle, _Anm);
}


std::tuple<ndarray,Cndarray>
PyBnm_integrand(int    n,
                int    m,
                int    sampling,
                double k,
                double w0,
                Vec    Offset,
                double Tolerance)
{

  argument args0;
           args0.n      = n;
           args0.m      = m;
           args0.rhon   = (n + 0.5);
           args0.Offset = Offset;
           args0.k      = k;
           args0.R0     = sqrt(Offset[0]*Offset[0] + Offset[1]*Offset[1]);
           args0.xi     = acos(Offset[0]/args0.R0);
           args0.s      = 1./(k*w0);
           args0.w0     = w0;
           args0.r      = args.rhon/k;

   Cndarray _Bnm = Cndarray(sampling);
   auto _Bnm_data = _Bnm.mutable_data();

   Cndarray Angle = Cndarray(sampling);
   auto Angle_data = Angle.mutable_data();

   Vec X = Linspace( 0.0, pi<double>(),sampling);

   for (auto i=0; i<sampling;i++)
   {
     _Bnm_data[i] = Bnm_integrand(X[i], args);
     Angle_data[i] = X[i];
   }
   return std::make_tuple(Angle, _Bnm);
 }





PYBIND11_MODULE(GaussianBeam, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

    module.def("Anm",
               &Anm,
               py::arg("n"),
               py::arg("m"),
               py::arg("k"),
               py::arg("w0"),
               py::arg("Offset"),
               py::arg("Tolerance") = 1e-8,
               "Compute Anm");


     module.def("Anm_integrand",
                &PyAnm_integrand,
                py::arg("n"),
                py::arg("m"),
                py::arg("sampling"),
                py::arg("k"),
                py::arg("w0"),
                py::arg("Offset"),
                py::arg("Tolerance") = 1e-8,
                "Compute Anm integrand as a function of theta");



    module.def("Bnm",
               &Bnm,
               py::arg("n"),
               py::arg("m"),
               py::arg("k"),
               py::arg("w0"),
               py::arg("Offset"),
               py::arg("Tolerance") = 1e-8,
               "Compute Bnm");


     module.def("Bnm_integrand",
                &PyBnm_integrand,
                py::arg("n"),
                py::arg("m"),
                py::arg("sampling"),
                py::arg("k"),
                py::arg("w0"),
                py::arg("Offset"),
                py::arg("Tolerance") = 1e-8,
                "Compute Bnm integrand as a function of theta");
}














//-
