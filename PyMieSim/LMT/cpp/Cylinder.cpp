//#include "Functions.cpp"
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#define j complex128(0.0,1.0)


namespace py = pybind11;

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;



namespace Cylinder
{



void
CoefficientAnBn(const double &Diameter,  // ref: https://doi.org/10.1364/AO.44.002338
                const double &Wavelength,
                const double &Index,
                const double &nMedium,
                const int    &MaxOrder,
                complex128   *an,
                complex128   *bn)
{

  double x = PI * Diameter /Wavelength;

  double mt = Index, m = nMedium;

  complex128 numerator, denominator;

  for (auto order = 1; order < MaxOrder+1; order++)
  {

    numerator   = mt * Jn(order, mt*x) * Jn_p(order, m*x) - m * Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = mt * Jn(order, mt*x) * Hn_p(order, m*x) - m * Jn_p(order, mt*x) * Hn(order, m*x);

    an[order] = numerator/denominator;

    numerator   = m * Jn(order, mt*x) * Jn_p(order, m*x) - mt*Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = m * Jn(order, mt*x) * Hn_p(order, m*x) - mt*Jn_p(order, mt*x) * Hn(order, m*x);
    bn[order] = numerator/denominator;
  }
}




void
MiePiTau(const double  mu,
         const int     OrderMax,
         complex128*   pin,
         complex128*   taun )

{
  pin[0] = 1.;
  pin[1] = 3. * mu;

  taun[0] = mu;
  taun[1] = 3.0 * cos(2. * acos(mu) );

  double n = 0;
  for (long unsigned int i = 2; i < OrderMax; i++)
      {
       n = (double)i;

       pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

       taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
     }
}


static double
C_Qsca(complex128    *an,
       complex128    *bn,
       const double   SizeParam,
       const Vec      n)

{
     int MaxOrder = GetMaxOrder(SizeParam);
     double Qsca = 2. / (SizeParam * SizeParam);
     complex128 temp;
     for(auto it = 0; it < MaxOrder; ++it)
     {
       temp += (2.*n[it]+1) * (   std::real( an[it] ) * std::real( an[it] )
                                + std::imag( an[it] ) * std::imag( an[it] )
                                + std::real( bn[it] ) * std::real( bn[it] )
                                + std::imag( bn[it] ) * std::imag( bn[it] ) );
     }
     return Qsca * std::real(temp);
}

static double
C_Qext(complex128    *an,
       complex128    *bn,
       const double   SizeParam,
       const Vec      n)

{
     int MaxOrder = GetMaxOrder(SizeParam);
     double Qsca = 2. / (SizeParam * SizeParam);
     complex128 temp;
     for(auto it = 0; it < MaxOrder; ++it)
     {
       temp += (2.*n[it]+1) * ( std::real( an[it] + an[it] ) );
     }
     return Qsca * std::real(temp);
}


std::tuple<double, double, double>
Efficiencies(const double  Diameter,
             const double  Wavelength,
             const double  Index,
             const double  nMedium)
{

    double Qsca,
           Qext, 
           Qabs,
           SizeParam = PI * Diameter /Wavelength;

    Vec n, n2;

    int MaxOrder = GetMaxOrder(PI * Diameter / Wavelength);

    std::tie(n, n2) = Arrange(1, MaxOrder + 1);

    complex128 *an     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn     = (complex128*) malloc(sizeof(complex128)*MaxOrder);

    CoefficientAnBn(Diameter, Wavelength, Index, nMedium, MaxOrder, an, bn);

    Qsca = C_Qsca(an, bn, SizeParam, n);

    Qext = C_Qext(an, bn, SizeParam, n);

    Qabs = Qext - Qsca;
    free(an);
    free(bn);
    return std::make_tuple(Qsca, Qext, Qabs);
}


static void
GetS1S2(const double          Index,
      const double            Diameter,
      const double            Wavelength,
      const double            nMedium,
      double                 *PhiPtr,
      const long unsigned int lenght,
      complex128*             s1s2)

{
    Vec n, n2;

    int MaxOrder = GetMaxOrder(PI * Diameter / Wavelength);

    std::tie(n, n2) = Arrange(1, MaxOrder + 1);

    complex128 *an     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *pin    = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *taun   = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *temp0  = &s1s2[0],
               *temp1  = &s1s2[lenght] ;

    CoefficientAnBn(Diameter, Wavelength, Index, nMedium, MaxOrder, an, bn);

    iVec S1   = iVec(lenght),
         S2   = iVec(lenght);

    for (long unsigned int i = 0; i < lenght; i++){

        MiePiTau(cos( PhiPtr[i]-PI/2 ), MaxOrder, pin, taun);

        for (auto k = 0; k < MaxOrder ; k++){
            *temp0 += n2[k] * ( an[k] * pin[k] +  bn[k] * taun[k] );
            *temp1 += n2[k] * ( an[k] * taun[k] + bn[k] * pin[k] ) ;

          }
    temp0++ ;
    temp1++ ;
    }

free(an);
free(bn);
free(pin);
free(taun);
}



std::tuple<Cndarray, Cndarray>
S1S2(const double            Index,
     const double            Diameter,
     const double            Wavelength,
     const double            nMedium,
     ndarray                 Phi)

{
    info PhiBuffer   = Phi.request();

    double *PhiPtr   = (double *) PhiBuffer.ptr;

    int lenght   = PhiBuffer.size,
        MaxOrder = GetMaxOrder(PI * Diameter / Wavelength);

    Vec n, n2;

    complex128 *an     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *pin    = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *taun   = (complex128*) malloc(sizeof(complex128)*MaxOrder),
                temp0  = 0.,
                temp1   = 0.;

    Cndarray S1 = ndarray(lenght,0),
             S2 = ndarray(lenght,0);

    auto S1_data = S1.mutable_data(),
         S2_data = S2.mutable_data();

    std::tie(n, n2) = Arrange(1, MaxOrder + 1);

    CoefficientAnBn(Diameter, Wavelength, Index, nMedium, MaxOrder, an, bn);

    for (auto i = 0; i < lenght; i++){

        MiePiTau(cos( PhiPtr[i]-PI/2 ), MaxOrder, pin, taun);

        for (auto k = 0; k < MaxOrder ; k++){
            temp0 += n2[k] * ( an[k] * pin[k] +  bn[k] * taun[k] );
            temp1 += n2[k] * ( an[k] * taun[k] + bn[k] * pin[k]  );
          }
        S1_data[i] = temp0;
        S2_data[i] = temp1;

        temp0 = 0.; temp1=0.;

    }
    free(an);
    free(bn);
    free(pin);
    free(taun);
    return std::make_tuple(S1,S2);
}





std::pair<Cndarray, Cndarray>
FieldsStructured(double     Index,
                   double     Diameter,
                   double     wavelength,
                   double     nMedium,
                   ndarray    Phi,
                   ndarray    Theta,
                   double     Polarization,
                   double     E0,
                   double     R)
{

  info PhiBuffer   = Phi.request();
  info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0],
      w            = 0;

  Cndarray result0 = Cndarray(PhiLength * ThetaLength);
  Cndarray result1 = Cndarray(PhiLength * ThetaLength);

  info buf0 = result0.request();
  info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          k = 2. * PI / wavelength,
          temp0 ;

  complex128 *EPhiPtr    = (complex128 *) buf0.ptr,
             *EThetaPtr  = (complex128 *) buf1.ptr,
              propagator  = E0 / (k * R) * exp(-J*k*R),
             *s1s2       = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopStructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, true, Polarization);


  result0.resize({PhiLength,ThetaLength});
  result1.resize({PhiLength,ThetaLength});

  free(s1s2);

  return std::make_pair(result0,result1);

}




std::pair<Cndarray, Cndarray>
FieldsUnstructured(double     Index,
                   double     Diameter,
                   double     wavelength,
                   double     nMedium,
                   ndarray    Phi,
                   ndarray    Theta,
                   double     Polarization,
                   double     E0,
                   double     R)
{

  info PhiBuffer   = Phi.request();
  info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0];

  Cndarray result0 = Cndarray(Phi.size());
  Cndarray result1 = Cndarray(Phi.size());

  info buf0 = result0.request();
  info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
         k         = 2. * PI / wavelength,
         temp0     = 0.;

  complex128 *EPhiPtr    = (complex128 *) buf0.ptr,
             *EThetaPtr  = (complex128 *) buf1.ptr,
              propagator = E0 / (k * R) * exp(-J*k*R),
             *s1s2       = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopUnstructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, true, Polarization);

  free(s1s2);
  free(EPhiPtr);
  free(EThetaPtr);

  return std::make_pair(result0.attr("transpose")(), result1.attr("transpose")());
}


std::pair<Cndarray, Cndarray>
FieldsStructuredUnpolarized(double     Index,
                            double     Diameter,
                            double     wavelength,
                            double     nMedium,
                            ndarray    Phi,
                            ndarray    Theta,
                            double     E0,
                            double     R)
{

  info PhiBuffer   = Phi.request();
  info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0],
      w            = 0;

  Cndarray result0 = Cndarray(PhiLength * ThetaLength);
  Cndarray result1 = Cndarray(PhiLength * ThetaLength);

  info buf0 = result0.request();
  info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          k        = 2. * PI / wavelength,
          temp0    = 0.,
          temp1    = 1./sqrt(2.) ;

  complex128 *EPhiPtr    = (complex128 *) buf0.ptr,
             *EThetaPtr  = (complex128 *) buf1.ptr,
              propagator = E0 / (k * R) * exp(-J*k*R),
             *s1s2       = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopStructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, false, 0.);

  free(s1s2);
  result0.resize({PhiLength,ThetaLength});
  result1.resize({PhiLength,ThetaLength});

  return std::make_pair(result0.attr("transpose")(), result1.attr("transpose")());
}





std::pair<Cndarray, Cndarray>
FieldsUnstructuredUnpolarized(double     Index,
                              double     Diameter,
                              double     wavelength,
                              double     nMedium,
                              ndarray    Phi,
                              ndarray    Theta,
                              double     E0,
                              double     R)
{

  info PhiBuffer   = Phi.request();
  info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0];

  Cndarray result0 = Cndarray(Phi.size());
  Cndarray result1 = Cndarray(Phi.size());

  info buf0 = result0.request();
  info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          k        = 2. * PI / wavelength,
          temp0    = 0.,
          temp1    = 1./sqrt(2.) ;

  complex128 *EPhiPtr    = (complex128 *) buf0.ptr,
             *EThetaPtr  = (complex128 *) buf1.ptr,
              propagator = E0 / (k * R) * exp(-J*k*R),
             *s1s2       = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopUnstructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, false, 0.);

  free(s1s2);

  return std::make_pair(result0.attr("transpose")(), result1.attr("transpose")());
}





}




// -
