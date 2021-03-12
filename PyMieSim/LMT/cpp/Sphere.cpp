#include "utils.cpp"
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;
typedef py::buffer_info info;


namespace Sphere
{




void
LowFrequencyMie_ab(const double Index,
                   const double SizeParam,
                   complex128   *an,
                   complex128   *bn)
{
  double LL, m2, x3, x4, x5, x6;
  complex128 a1, a2, b1, b2;

  m2 = Index * Index;
  LL = (m2 - 1) / (m2 + 2);
  x3 = SizeParam * SizeParam * SizeParam;
  x4 = x3 * SizeParam;
  x5 = x4 * SizeParam;
  x6 = x5 * SizeParam;

  a1 = (-2.*J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  a2 = (-1.*J * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  b1 = (-1.*J * x5 / 45.) * (m2 - 1.);
  b2 = 0. + 0.*J;

  an[0] = a1;
  an[1] = a2;
  bn[0] = b1;
  bn[1] = b2;
  return ;
}



void
HighFrequencyMie_ab(const double               m,
                    const double               x,
                    const long unsigned int    OrderMax,
                    const                      Vec n,
                    complex128*                an,
                    complex128*                bn)

{
  const double mx = m * x;
  const double temp  = sqrt(0.5 * PI * x);
  const long unsigned int nmx = (long unsigned int) ( std::max( OrderMax, (long unsigned int) abs(mx) ) + 16 );
  iVec gsx, gs1x;
  iVec px, chx, p1x, ch1x, D, da, db;
  Vec Dn = Vec(nmx);

  p1x.push_back( sin(x) );
  ch1x.push_back( cos(x) );

  for (double i = nmx - 1; i > 1; i--)
  {
      Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );
  }

  for (long unsigned int i = 0; i < OrderMax; i++)
  {
    px.push_back(  temp * jn( n[i], x ) );
    chx.push_back(-temp * yn( n[i], x ) );

    p1x.push_back(px[i]);
    ch1x.push_back(chx[i]);

    gsx.push_back( px[i] - 1.*J * chx[i] );
    gs1x.push_back( p1x[i] - 1.*J * ch1x[i] );

    D.push_back(Dn[i+1]);

    da.push_back( D[i] / m + n[i] / x );
    db.push_back( m * D[i] + n[i] / x );

    an[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;
    bn[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;
  }
}




void
CoefficientAnBn(const double &Diameter,
                const double &Wavelength,
                const double &Index,
                const double &nMedium,
                const int    &MaxOrder,
                complex128   *an,
                complex128   *bn)
{
  Vec n, n2;
  double SizeParam = GetSizeParameter(Diameter, Wavelength, nMedium);

  std::tie(n, n2) = Arrange(1, MaxOrder + 1);

  if (SizeParam < 0.5){LowFrequencyMie_ab(Index, SizeParam, an, bn); }

  else{HighFrequencyMie_ab(Index, SizeParam, MaxOrder, n, an, bn);}
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
  for (auto i = 2; i < OrderMax; i++)
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

    Vec n, n2;

    double SizeParam = GetSizeParameter(Diameter, Wavelength, nMedium),
           Qsca      = 0.,
           Qext      = 0.,
           Qabs      = 0.;

    int MaxOrder = GetMaxOrder(SizeParam);

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
GetS1S2(const double            Index,
        const double            Diameter,
        const double            Wavelength,
        const double            nMedium,
        double                 *PhiPtr,
        const long unsigned int lenght,
        complex128*             s1s2)

{

    iVec S1   = iVec(lenght),
         S2   = iVec(lenght);

    int MaxOrder = GetMaxOrder( GetSizeParameter(Diameter, Wavelength, nMedium) );

    complex128 *an        = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn        = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *pin       = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *taun      = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *temp0     = &s1s2[0],
               *temp1     = &s1s2[lenght],
                prefactor = 0.;

    CoefficientAnBn(Diameter, Wavelength, Index, nMedium, MaxOrder, an, bn);

    for (long unsigned int i = 0; i < lenght; i++){

        MiePiTau(cos( PhiPtr[i]-PI/2 ), MaxOrder, pin, taun);

        for (auto k = 0; k < MaxOrder ; k++){
             prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            *temp0    += prefactor * ( an[k] * pin[k] +  bn[k] * taun[k] );
            *temp1    += prefactor * ( an[k] * taun[k] + bn[k] * pin[k] ) ;

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
S1S2(const double  Index,
     const double  Diameter,
     const double  Wavelength,
     const double  nMedium,
     ndarray       Phi)

{
    info        PhiBuffer = Phi.request();

    double     *PhiPtr    = (double *) PhiBuffer.ptr;

    int         lenght    = PhiBuffer.size,
                MaxOrder  = GetMaxOrder( GetSizeParameter(Diameter, Wavelength, nMedium) );

    Vec         n,
                n2;

    complex128 *an      = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn      = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *pin     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *taun    = (complex128*) malloc(sizeof(complex128)*MaxOrder),
                temp0   = 0.,
                temp1   = 0.;

    Cndarray   S1       = ndarray(lenght,0),
               S2       = ndarray(lenght,0);

    auto       S1_data  = S1.mutable_data(),
               S2_data  = S2.mutable_data();

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
                   double     Wavelength,
                   double     nMedium,
                   ndarray    Phi,
                   ndarray    Theta,
                   double     Polarization,
                   double     E0,
                   double     R)
{



  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0],
              w            = 0;

  Cndarray    EPhi         = Cndarray(PhiLength * ThetaLength),
              ETheta       = Cndarray(PhiLength * ThetaLength);

  info        buf0         = EPhi.request(),
              buf1         = ETheta.request();

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr,
              k            = 2. * PI / Wavelength;

  complex128 *EPhiPtr      = (complex128 *) buf0.ptr,
             *EThetaPtr    = (complex128 *) buf1.ptr,
              propagator   = E0 / (k * R) * exp(-J*k*R),
             *s1s2         = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, Wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopStructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, true, Polarization);

  EPhi.resize({PhiLength,ThetaLength}); ETheta.resize({PhiLength,ThetaLength});

  free(s1s2);

  return std::make_pair(EPhi,ETheta);

}




std::pair<Cndarray, Cndarray>
FieldsUnstructured(double     Index,
                   double     Diameter,
                   double     Wavelength,
                   double     nMedium,
                   ndarray    Phi,
                   ndarray    Theta,
                   double     Polarization,
                   double     E0,
                   double     R)
{

  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0],
              w            = 0;

  Cndarray    EPhi         = Cndarray(PhiLength * ThetaLength),
              ETheta       = Cndarray(PhiLength * ThetaLength);

  info        buf0         = EPhi.request(),
              buf1         = ETheta.request();

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr,
              k            = 2. * PI / Wavelength;

  complex128 *EPhiPtr      = (complex128 *) buf0.ptr,
             *EThetaPtr    = (complex128 *) buf1.ptr,
              propagator   = E0 / (k * R) * exp(-J*k*R),
             *s1s2         = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, Wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopUnstructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, true, Polarization);

  free(s1s2);

  return std::make_pair(EPhi.attr("transpose")(), ETheta.attr("transpose")());
}


std::pair<Cndarray, Cndarray>
FieldsStructuredUnpolarized(double     Index,
                            double     Diameter,
                            double     Wavelength,
                            double     nMedium,
                            ndarray    Phi,
                            ndarray    Theta,
                            double     E0,
                            double     R)
{

  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0],
              w            = 0;

  Cndarray    EPhi         = Cndarray(PhiLength * ThetaLength),
              ETheta       = Cndarray(PhiLength * ThetaLength);

  info        buf0         = EPhi.request(),
              buf1         = ETheta.request();

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr,
              k            = 2. * PI / Wavelength;

  complex128 *EPhiPtr      = (complex128 *) buf0.ptr,
             *EThetaPtr    = (complex128 *) buf1.ptr,
              propagator   = E0 / (k * R) * exp(-J*k*R),
             *s1s2         = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, Wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopStructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, false, 0.);

  EPhi.resize({PhiLength,ThetaLength}); ETheta.resize({PhiLength,ThetaLength});

  free(s1s2);

  return std::make_pair(EPhi.attr("transpose")(), ETheta.attr("transpose")());
}





std::pair<Cndarray, Cndarray>
FieldsUnstructuredUnpolarized(double     Index,
                              double     Diameter,
                              double     Wavelength,
                              double     nMedium,
                              ndarray    Phi,
                              ndarray    Theta,
                              double     E0,
                              double     R)
{

  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0];

  Cndarray    EPhi         = Cndarray(PhiLength * ThetaLength),
              ETheta       = Cndarray(PhiLength * ThetaLength);

  info        buf0         = EPhi.request(),
              buf1         = ETheta.request();

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr,
              k            = 2. * PI / Wavelength;

  complex128 *EPhiPtr      = (complex128 *) buf0.ptr,
             *EThetaPtr    = (complex128 *) buf1.ptr,
              propagator   = E0 / (k * R) * exp(-J*k*R),
             *s1s2         = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  GetS1S2(Index, Diameter, Wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  LoopUnstructured(PhiLength, ThetaLength, propagator, PhiPtr, ThetaPtr, EPhiPtr, EThetaPtr, s1s2, false, 0.);

  free(s1s2);

  return std::make_pair(EPhi.attr("transpose")(), ETheta.attr("transpose")());
}




}




// -
