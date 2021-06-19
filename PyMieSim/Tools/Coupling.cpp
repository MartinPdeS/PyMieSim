#include <vector>
#include <complex>
#include <cmath>
#include <ostream>
#include <tuple>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

typedef std::complex<double> complex128;
typedef py::array_t<double> ndarray;
typedef py::array_t<std::complex<double>> Cndarray;
typedef py::buffer_info info;
using namespace std;


std::tuple<complex128, complex128>
CoherentLoop(Cndarray& ScalarField, Cndarray& ETheta, Cndarray& EPhi)
{
  complex128   CouplingTheta  = 0.0,
               CouplingPhi    = 0.0;

  complex128 * ScalarFieldPtr = (complex128 *) ScalarField.request().ptr,
             * EPhiPtr        = (complex128 *) EPhi.request().ptr,
             * EThetaPtr      = (complex128 *) ETheta.request().ptr;

  uint size = ScalarField.request().size ;

  for (uint i=0; i<size; i++)
  {
    CouplingTheta += ScalarFieldPtr[i] * EThetaPtr[i];
    CouplingPhi   += ScalarFieldPtr[i] * EPhiPtr[i];
  }

  return tie(CouplingTheta, CouplingPhi);
}

std::tuple<double, double>
NoCoherentLoop(Cndarray& ScalarField, Cndarray& ETheta, Cndarray& EPhi)
{
  double       CouplingTheta  = 0.0,
               CouplingPhi    = 0.0;

  complex128 * ScalarFieldPtr = (complex128 *) ScalarField.request().ptr,
             * EPhiPtr        = (complex128 *) EPhi.request().ptr,
             * EThetaPtr      = (complex128 *) ETheta.request().ptr;

  uint size = ScalarField.request().size ;

  for (uint i=0; i<size; i++)
  {
    CouplingTheta += pow( abs( ScalarFieldPtr[i] * EThetaPtr[i]), 2 );
    CouplingPhi   += pow( abs( ScalarFieldPtr[i] * EPhiPtr[i]),   2 );
  }

  return tie(CouplingTheta, CouplingPhi);
}

//_________________________POINT_COUPLING_________________________________________________________________________________

double
NoCoherentPointCouplingFilter(Cndarray&    ScalarField,
                              Cndarray&    ETheta,
                              Cndarray&    EPhi,
                              const double dOmega,
                              const double Filter)
{
  double CouplingTheta,
         CouplingPhi,
         ThetaFiltering  = pow( sin(Filter), 2 ),
         PhiFiltering    = pow( cos(Filter), 2 );

  std::tie(CouplingTheta, CouplingPhi) = NoCoherentLoop(ScalarField, ETheta, EPhi) ;

  CouplingTheta *= ThetaFiltering;
  CouplingPhi   *= PhiFiltering;

  return CouplingTheta + CouplingPhi;
}


double
NoCoherentPointCoupling(Cndarray&    ScalarField,
                        Cndarray&    ETheta,
                        Cndarray&    EPhi,
                        const double dOmega)
{
  double CouplingTheta, CouplingPhi;

  std::tie(CouplingTheta, CouplingPhi) = NoCoherentLoop(ScalarField, ETheta, EPhi) ;

  return (CouplingTheta + CouplingPhi) * dOmega;
}




double
CoherentPointCouplingFilter(Cndarray&    ScalarField,
                            Cndarray&    ETheta,
                            Cndarray&    EPhi,
                            const double dOmega,
                            const double Filter)
{
  complex128   CouplingTheta, CouplingPhi;

  double       ThetaFiltering = pow( sin(Filter), 2 ),
               PhiFiltering = pow( cos(Filter), 2 );

  std::tie(CouplingTheta, CouplingPhi) = CoherentLoop(ScalarField, ETheta, EPhi) ;

  CouplingTheta = pow( abs(CouplingTheta), 2 )   * ThetaFiltering;
  CouplingPhi   = pow( abs(CouplingPhi),   2 )   * PhiFiltering;

  return  abs(CouplingTheta + CouplingPhi) * dOmega;
}


double
CoherentPointCoupling(Cndarray&    ScalarField,
                      Cndarray&    ETheta,
                      Cndarray&    EPhi,
                      const double dOmega)
{
  complex128 CouplingTheta = 0.0, CouplingPhi   = 0.0;

  std::tie(CouplingTheta, CouplingPhi) = CoherentLoop(ScalarField, ETheta, EPhi) ;

  CouplingTheta = pow( abs(CouplingTheta), 2 );
  CouplingPhi   = pow( abs(CouplingPhi), 2 );

  return abs(CouplingTheta + CouplingPhi) * dOmega;
}




//_________________________Mean_COUPLING_________________________________________________________________________________

double
CoherentMeanCoupling(Cndarray&     ScalarField,
                     Cndarray&     ETheta,
                     Cndarray&     EPhi,
                     const double  dOmega,
                     const double  Omega)
{
  double CouplingTheta, CouplingPhi;

  double factor = dOmega / Omega;

  std::tie(CouplingTheta, CouplingPhi) = NoCoherentLoop(ScalarField, ETheta, EPhi) ;

  CouplingTheta = CouplingTheta;
  CouplingPhi   = CouplingPhi;

  return abs( CouplingTheta + CouplingPhi ) * factor;
}


double
CoherentMeanCouplingFilter(Cndarray&     ScalarField,
                           Cndarray&     ETheta,
                           Cndarray&     EPhi,
                           const double  dOmega,
                           const double  Omega,
                           const double  Filter)
{
  double CouplingTheta,
         CouplingPhi;

  double ThetaFiltering = pow( sin(Filter), 2 ),
         PhiFiltering   = pow( cos(Filter), 2 );

  std::tie(CouplingTheta, CouplingPhi) = NoCoherentLoop(ScalarField, ETheta, EPhi) ;

  CouplingTheta = CouplingTheta * ThetaFiltering;
  CouplingPhi   = CouplingPhi   * PhiFiltering;

  return abs( CouplingTheta + CouplingPhi ) * dOmega / Omega;
}


double
NoCoherentMeanCoupling(Cndarray&    ScalarField,
                       Cndarray&    ETheta,
                       Cndarray&    EPhi,
                       const double dOmega,
                       const double Omega)
{
  return NoCoherentPointCoupling(ScalarField, ETheta,  EPhi, dOmega);
}


double
NoCoherentMeanCouplingFilter(Cndarray&     ScalarField,
                             Cndarray&     ETheta,
                             Cndarray&     EPhi,
                             const double  dOmega,
                             const double  Omega,
                             const double  Filter)
{
  return NoCoherentPointCouplingFilter(ScalarField, ETheta, EPhi, dOmega, Filter);
}



PYBIND11_MODULE(_Coupling, module) {
    module.doc() = "Coherent and non-coherent coupling";

    module.def("NoCoherentPointCouplingFilter", &NoCoherentPointCouplingFilter,
               py::arg("ScalarField"),
               py::arg("ETheta"),
               py::arg("EPhi"),
               py::arg("dOmega"),
               py::arg("Filter") );

    module.def("NoCoherentPointCoupling", &NoCoherentPointCoupling,
              py::arg("ScalarField"),
              py::arg("ETheta"),
              py::arg("EPhi"),
              py::arg("dOmega"));


    module.def("NoCoherentMeanCouplingFilter", &NoCoherentMeanCouplingFilter,
               py::arg("ScalarField"),
               py::arg("ETheta"),
               py::arg("EPhi"),
               py::arg("dOmega"),
               py::arg("Omega"),
               py::arg("Filter") );

    module.def("NoCoherentMeanCoupling", &NoCoherentMeanCoupling,
              py::arg("ScalarField"),
              py::arg("ETheta"),
              py::arg("EPhi"),
              py::arg("dOmega"),
              py::arg("Omega") );


    module.def("CoherentPointCouplingFilter", &CoherentPointCouplingFilter,
               py::arg("ScalarField"),
               py::arg("ETheta"),
               py::arg("EPhi"),
               py::arg("dOmega"),
               py::arg("Filter") );

    module.def("CoherentPointCoupling", &CoherentPointCoupling,
              py::arg("ScalarField"),
              py::arg("ETheta"),
              py::arg("EPhi"),
              py::arg("dOmega"));


    module.def("CoherentMeanCouplingFilter", &CoherentMeanCouplingFilter,
               py::arg("ScalarField"),
               py::arg("ETheta"),
               py::arg("EPhi"),
               py::arg("dOmega"),
               py::arg("Omega"),
               py::arg("Filter") );

    module.def("CoherentMeanCoupling", &CoherentMeanCoupling,
              py::arg("ScalarField"),
              py::arg("ETheta"),
              py::arg("EPhi"),
              py::arg("dOmega"),
              py::arg("Omega") );

}
























// -
