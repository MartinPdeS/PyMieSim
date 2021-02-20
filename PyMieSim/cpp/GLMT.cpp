#include <boost/math/special_functions/bessel_prime.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#include "Math.cpp"
#include "Special.cpp"
#include <cmath>
#include <tuple>

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<std::vector<double>> Matrix;




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
  printf("\n\n\n %f\n \n\n\n", (jnm[0] * yn[1] - mt * yn[0] * jnm[1]) );
  complex128 b = (jnm[0] * jn[1] - mt * jn[0] * jnm[1]) / (jnm[0] * yn[1] - mt * yn[0] * jnm[1]);



  for (long unsigned int i = 0; i < Order; i++){
    an.push_back(a);
    bn.push_back(b);
  }
  //printVector(an);
  //printVector(bn);

  return std::make_tuple(an, bn); 
}









std::tuple<iVec*, iVec*, iVec*, Vec, Vec, Vec>
ScatteredField(double      Radius,
               std::size_t Order,
               double      Eps,
               double      Mu,
               double      Wavelength,
               int         Num)
{

  double k = 2 * PI / Wavelength;
  iVec an, bn;
  complex128 j (0., 1.0);
  iVec *RCompE,
       *RCompO,
       *PhiCompE,
       *PhiCompO,
       *ThetaCompE,
       *ThetaCompO;

  Vec R, Phi, Theta;


  std::tie(an, bn) = SphereCoefficient(Radius,
                                       Order,
                                       Eps,
                                       Mu,
                                       Wavelength);

  int lmax = an.size();


  //std::tie(R, Phi, Theta) = Fibonacci_sphere(10, 0.3);
  std::tie(R, Phi, Theta) = FullMesh(Num);



  iVec* E0 = new iVec(R.size());
  iVec* E1 = new iVec(R.size());
  iVec* E2 = new iVec(R.size());



  for (int L = 1; L < lmax+1; L++){


      complex128 L_i = (complex128) L;
      complex128 En = pow(j,L_i) * (2.*L_i+1.) / ( L_i*(L_i+1.) );


      std::tie(RCompE, PhiCompE, ThetaCompE) = N3e1n(L,
                                                     k,
                                                     R,
                                                     Phi,
                                                     Theta);


     std::tie(RCompO, PhiCompO, ThetaCompO) = M3o1n(L,
                                                    k,
                                                    R,
                                                    Phi,
                                                    Theta);


      for (long unsigned int l = 0; l < R.size(); l++){
        (*E0)[l] += En*(j * an[L-1] * (*RCompE)[l]     - bn[L-1] * (*RCompO)[l] );
        (*E1)[l] += En*(j * an[L-1] * (*PhiCompE)[l]   - bn[L-1] * (*PhiCompO)[l] );
        (*E2)[l] += En*(j * an[L-1] * (*ThetaCompE)[l] - bn[L-1] * (*ThetaCompO)[l]);

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
