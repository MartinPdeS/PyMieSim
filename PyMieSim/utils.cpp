#include <iostream>
#include <vector>
#include <complex>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

typedef std::complex<double> complex128;
typedef std::vector<double> Vec;
typedef py::array_t<double> ndarray;

#define J complex128(0.0,1.0)
#define PI (double) 3.14159265358979323846264338

using namespace std;



void
Cart2sph(Vec     X,
         Vec     Y,
         Vec     Z,
         double* R,
         double* Phi,
         double* Theta){

  for (auto i = 0; i < X.size(); i++){
    R[i]     = sqrt( pow( X[i],2) + pow(Y[i],2) + pow(Z[i],2) ) ;
    Phi[i]   = atan( Y[i]/X[i] );
    Theta[i] = atan( sqrt( pow(Y[i],2) / pow(X[i],2) ) );
  }
  return ;
}


std::tuple<ndarray, ndarray>
Fibonacci_Spherical(int Samples, double MaxAngle)
{
    Vec X, Y, Z;

    double phi          = PI * (3. - sqrt(5.)),
           SolidAngle   = abs( 2*PI * (cos(MaxAngle) - 1)),
           ratio        = 4*PI / SolidAngle;

    int    TrueSampling = (int)( Samples * ratio );

    for (auto i = 0; i < TrueSampling; i++){

        double y      = 1 - (i / (double)(TrueSampling - 1)) * 2,
               radius = sqrt(1 - y * y),
               theta  = phi * i,
               x      = cos(theta) * radius,
               z      = sin(theta) * radius;

        X.push_back(x);
        Y.push_back(y);
        Z.push_back(z);

        if (i >= (int)Samples - 1) break;
      }

    ndarray R          = ndarray(X.size()),
            Phi        = ndarray(X.size()),
            Theta      = ndarray(X.size());

    auto    R_data     = R.mutable_data(),
            Phi_data   = Phi.mutable_data(),
            Theta_data = Theta.mutable_data();

    Cart2sph(X,Y,Z, R_data, Phi_data, Theta_data);


    return std::make_tuple(Phi, Theta);
}



std::tuple<ndarray, ndarray, ndarray>
Fibonacci_Cartesian(int Samples, double MaxAngle)
{
    ndarray X            = ndarray(Samples),
            Y            = ndarray(Samples),
            Z            = ndarray(Samples);

    auto    X_data       = X.mutable_data(),
            Y_data       = Y.mutable_data(),
            Z_data       = Z.mutable_data();

    double  phi          = PI * (3. - sqrt(5.)),
            SolidAngle   = abs( 2 * PI * ( cos(MaxAngle) - 1 ) ),
            ratio        = 4 * PI / SolidAngle;

    int     TrueSampling = (int)( Samples * ratio );

    for (auto i = 0; i < TrueSampling; i++){

        double y      = 1 - (i / (double)(TrueSampling - 1)) * 2,
               radius = sqrt(1 - y * y),
               theta  = phi * i,
               x      = cos(theta) * radius,
               z      = sin(theta) * radius;


        X_data[i] = x ;
        Y_data[i] = y ;
        Z_data[i] = z ;

        if (i >= (int)Samples - 1) break;
      }


    return std::make_tuple(X, Y, Z);
}


PYBIND11_MODULE(_utils, module) {
    module.doc() = "LGeneralized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

    module.def("Fibonacci_Cartesian",
               &Fibonacci_Cartesian,
               py::arg("Samples"),
               py::arg("MaxAngle"),
               "Return S1");


     module.def("Fibonacci_Spherical",
                &Fibonacci_Spherical,
                py::arg("Samples"),
                py::arg("MaxAngle"),
                "Return S1");



}

























// --
