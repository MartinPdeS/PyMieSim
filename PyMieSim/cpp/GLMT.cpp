#include <boost/math/special_functions/bessel_prime.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <boost/geometry/algorithms/transform.hpp>
namespace py = pybind11;
#include "Math.cpp"
#include "Special.cpp"
#include <cmath>
#include <tuple>
typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<std::vector<complex128>> iMatrix;


std::tuple<iVec, iVec>
SphereCoefficient(double Radius,
                  std::size_t Order,
                  double Eps,
                  double Mu,
                  double Wavelength)
{
  double k = 2 * PI / Wavelength;
  double SizeParam = k * Radius;
  double m = pow(Eps, 0.5);
  double mt = m / Mu;
  iVec an, bn;

  Vec jn  = Riccati1(Order, SizeParam);
  iVec yn = Riccati3(Order, SizeParam);
  Vec jnm = Riccati1(Order, m * SizeParam);

  complex128 a = (mt * jnm[0] * jn[1] - jn[0] * jnm[1]) / (mt * jnm[0] * yn[1] - yn[0] * jnm[1]);
  complex128 b = (jnm[0] * jn[1] - mt * jn[0] * jnm[1]) / (jnm[0] * yn[1] - mt * yn[0] * jnm[1]);


  for (long unsigned int i = 0; i < Order; i++){
    an.push_back(a);
    bn.push_back(b);
  }

  return std::make_tuple(an, bn);
}



std::tuple<iVec*, iVec*, iVec*>
_ScatteredField(py::list an, py::list bn, double Wavelength, py::list l){
  double k = 2 * PI / Wavelength;
  int lmax = an.size();


  complex128 j (0., 1.0);
  iVec *RComp, *PhiComp, *ThetaComp;
  Vec R, Phi, Theta;
  R.push_back(1.);Phi.push_back(0.); Theta.push_back(0.);
  iVec* E0 = new iVec(R.size());
  iVec* E1 = new iVec(R.size());
  iVec* E2 = new iVec(R.size());


  for (int L = 1; L < lmax+1; L++){

      VectorSphericalHarmonics3 VSH = VectorSphericalHarmonics3(L, Wavelength);
      complex128 L_i = (complex128) L;
      complex128 En = pow(j,L_i) * (2.*L_i+1.) / ( L_i*(L_i+1.) );


      std::tie(RComp, PhiComp, ThetaComp) = VSH.N_e1n(R, Phi, Theta);

      //complex128 a = an.cast<float>();
      //printf("%f\n", a);

      for (long unsigned int l = 0; l < R.size(); l++){

        //(*E0)[l] = ((double)an[0]);// En*(j * an[L-1] * (*RComp)[l]     - bn[L-1] * (*RComp)[l] );
        //(*E1)[l] += En*(j * an[L-1] * (*PhiComp)[l]   - bn[L-1] * (*PhiComp)[l] );
        //(*E2)[l] += En*(j * an[L-1] * (*ThetaComp)[l] - bn[L-1] * (*ThetaComp)[l]);

    }
  }


  return std::make_tuple(E0, E1, E2);
}








std::tuple<iVec*, iVec*, iVec*, Vec, Vec, Vec>
ScatteredField(double      Radius,
               std::size_t Order,
               double      Eps,
               double      Mu,
               double      Wavelength,
               py::list    TTheta){

  double k = 2 * PI / Wavelength;
  iVec an, bn;
  complex128 j (0., 1.0);




  std::tie(an, bn) = SphereCoefficient(Radius,
                                       Order,
                                       Eps,
                                       Mu,
                                       Wavelength);

  int lmax = an.size();


  iVec *RComp, *PhiComp, *ThetaComp;
  Vec R, Phi, Theta;


  std::tie(R, Phi, Theta) = Fibonacci_sphere(100, 0.3);


  iVec* E0 = new iVec(R.size());
  iVec* E1 = new iVec(R.size());
  iVec* E2 = new iVec(R.size());



  for (int L = 1; L < lmax+1; L++){

      VectorSphericalHarmonics3 VSH = VectorSphericalHarmonics3(L, Wavelength);
      complex128 L_i = (complex128) L;
      complex128 En = pow(j,L_i) * (2.*L_i+1.) / ( L_i*(L_i+1.) );



      std::tie(RComp, PhiComp, ThetaComp) = VSH.N_e1n(R, Phi, Theta);


      for (long unsigned int l = 0; l < R.size(); l++){

        (*E0)[l] +=  En*(j * an[L-1] * (*RComp)[l]     - bn[L-1] * (*RComp)[l] );
        (*E1)[l] += En*(j * an[L-1] * (*PhiComp)[l]   - bn[L-1] * (*PhiComp)[l] );
        (*E2)[l] += En*(j * an[L-1] * (*ThetaComp)[l] - bn[L-1] * (*ThetaComp)[l]);

    }
  }


  return std::make_tuple(E0, E1, E2, R, Phi, Theta);
  }





PYBIND11_MODULE(GLMT, module) {
    py::module Sphere = module.def_submodule("Sphere", "Sphere scattering object");
    Sphere.def("Coefficient", &SphereCoefficient, "Return coefficient an & bn");
    Sphere.def("ScatteredField", &ScatteredField, "Return scattered field");

}
















//-
