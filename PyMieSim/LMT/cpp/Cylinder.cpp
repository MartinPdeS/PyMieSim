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





std::tuple<iVec, iVec>
LowFrequencyMie_ab(const double m,
                   const double x)
{
  iVec an, bn ;

  double LL, m2, x3, x4, x5, x6;
  complex128 a1, a2, b1, b2;

  m2 = m * m;
  LL = (m2 - 1) / (m2 + 2);
  x3 = x * x * x;
  x4 = x3 * x;
  x5 = x4 * x;
  x6 = x5 * x;

  a1 = (-2.*j * x3 / 3.) * LL - (2.*j * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  a2 = (-1.*j * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  b1 = (-1.*j * x5 / 45.) * (m2 - 1.);
  b2 = 0. + 0.*j;

  an.push_back(a1);
  an.push_back(a2);
  bn.push_back(b1);
  bn.push_back(b2);
  return std::make_tuple(an, bn);
}


void
CoefficientAnBn(const double &SizeParam,  // ref: https://doi.org/10.1364/AO.44.002338
             const double &Index,
             const double &nMedium,
             const int    &MaxOrder,
             complex128   *an,
             complex128   *bn)
{
  double alpha = SizeParam,
         beta  = alpha * Index,
         MuSp  = 1.,
         Mu    = 1.,
         M     = Index/nMedium;

  double mt = Index, m = nMedium, x = SizeParam;

  complex128 numerator, denominator, PsiAlpha, PsiBeta, PsiPBeta, PsiPAlpha, XiAlpha, XiPAlpha;

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



std::tuple<iVec, iVec>
MiePiTau(const double            mu,
         const long unsigned int MaxOrder)

{
  iVec pin = iVec(MaxOrder), taun = iVec(MaxOrder);
  pin[0] = 1.;
  pin[1] = 3. * mu;

  taun[0] = mu;
  taun[1] = 3.0 * cos(2. * acos(mu) );
  double n = 0;
  for (long unsigned int i = 2; i < MaxOrder; i++)
      {
       n = (double)i;

       pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

       taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
     }
     return std::make_tuple(pin, taun);
}




void
MiePiTau(const double            mu,
         const long unsigned int MaxOrder,
         iVec                    pin,
         iVec                    taun )

{
  pin[0] = 1.;
  pin[1] = 3. * mu;

  taun[0] = mu;
  taun[1] = 3.0 * cos(2. * acos(mu) );
  double n = 0;
  for (long unsigned int i = 2; i < MaxOrder; i++)
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
     for(long unsigned int it = 0; it < MaxOrder; ++it)
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
     for(long unsigned int it = 0; it < MaxOrder; ++it)
     {
       temp += (2.*n[it]+1) * ( std::real( an[it] + an[it] ) );
     }
     return Qsca * std::real(temp);
}



static void
_S1S2(const double            Index,
      const double            diameter,
      const double            wavelength,
      const double            nMedium,
      double                 *PhiPtr,
      const long unsigned int lenght,
      complex128*             s1s2)

{
    Vec n, n2;

    double m         = Index / nMedium,
           w         = wavelength / nMedium,
           SizeParam = PI * diameter / w;

    int MaxOrder = GetMaxOrder(SizeParam);

    std::tie(n, n2) = Arrange(1, MaxOrder + 1);

    complex128 *an     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn     = (complex128*) malloc(sizeof(complex128)*MaxOrder);

    CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

    iVec S1 = iVec(lenght),
         S2 = iVec(lenght),
         pin =iVec(MaxOrder), taun = iVec(MaxOrder);

    complex128 *temp0 = &s1s2[0],
               *temp1 = &s1s2[lenght] ;

    for (long unsigned int i = 0; i < lenght; i++){

        std::tie(pin, taun) = MiePiTau(cos( PhiPtr[i]-PI/2 ), MaxOrder);

        for (auto k = 0; k < MaxOrder ; k++){
            *temp0 += n2[k] * ( an[k] * pin[k] +  bn[k] * taun[k] );
            *temp1 += n2[k] * ( an[k] * taun[k] + bn[k] * pin[k] ) ;

          }
    temp0++ ;
    temp1++ ;
    }
    free(an);
    free(bn);

}






std::tuple<double, double, double>
Efficiencies(const double  Index, const double  SizeParam)

{

    double Qsca, Qext, Qabs, nMedium=1.;

    Vec n, n2;

    int MaxOrder = GetMaxOrder(SizeParam);

    std::tie(n, n2) = Arrange(1, MaxOrder + 1);

    complex128 *an     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn     = (complex128*) malloc(sizeof(complex128)*MaxOrder);

    CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

    Qsca = C_Qsca(an, bn, SizeParam, n);

    Qext = C_Qext(an, bn, SizeParam, n);

    Qabs = Qext - Qsca;
    free(an);
    free(bn);
    return std::make_tuple(Qsca, Qext, Qabs);
}


std::tuple<Cndarray, Cndarray>
S1S2(const double            Index,
     const double            Diameter,
     const double            wavelength,
     const double            nMedium,
     ndarray                 Phi)

{
    py::buffer_info PhiBuffer   = Phi.request();

    double *PhiPtr    = (double *) PhiBuffer.ptr,
            m         = Index / nMedium,
            w         = wavelength / nMedium,
            SizeParam = PI * Diameter / w;

    int lenght   = PhiBuffer.size,
        MaxOrder = GetMaxOrder(SizeParam);

    iVec pin  = iVec(MaxOrder),
         taun = iVec(MaxOrder);

    Vec n, n2;

    complex128 temp0=0., temp1=0.;

    Cndarray S1 = ndarray(lenght,0),
             S2 = ndarray(lenght,0);

    auto S1_data = S1.mutable_data(),
         S2_data = S2.mutable_data();

    std::tie(n, n2) = Arrange(1, MaxOrder + 1);

    complex128 *an     = (complex128*) malloc(sizeof(complex128)*MaxOrder),
               *bn     = (complex128*) malloc(sizeof(complex128)*MaxOrder);

    CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

    for (auto i = 0; i < lenght; i++){

        std::tie(pin, taun) = MiePiTau(cos( PhiPtr[i]-PI/2. ), MaxOrder);

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
    return std::make_tuple(S1,S2);
}



std::pair<Cndarray, Cndarray>
FieldsStructured(double     index,
                   double     diameter,
                   double     wavelength,
                   double     nMedium,
                   ndarray    Phi,
                   ndarray    Theta,
                   double     Polarization,
                   double     E0,
                   double     R)
{

  py::buffer_info PhiBuffer   = Phi.request();
  py::buffer_info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0],
      w            = 0;

  py::array_t<complex128> result0 = py::array_t<complex128>(PhiLength * ThetaLength);
  py::array_t<complex128> result1 = py::array_t<complex128>(PhiLength * ThetaLength);

  py::buffer_info buf0 = result0.request();
  py::buffer_info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          k = 2. * PI / wavelength,
          temp0 ;

  complex128 *ptr0 = (complex128 *) buf0.ptr,
             *ptr1 = (complex128 *) buf1.ptr,
             temp2,
             propagator = E0 / (k * R) * exp(-j*k*R),
             *s1s2 = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  _S1S2(index, diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          temp0 = ThetaPtr[t] ;
          ptr0[w] = j* propagator * s1s2[p] * (complex128) abs(cos(temp0 + Polarization));
          ptr1[w] =- propagator * s1s2[p + PhiLength] * (complex128) abs(sin(temp0 + Polarization));
          w++;
       }
  }

  free(s1s2);
  result0.resize({PhiLength,ThetaLength});
  result1.resize({PhiLength,ThetaLength});

  return std::make_pair(result0.attr("transpose")(), result1.attr("transpose")());
}






std::pair<Cndarray, Cndarray>
FieldsUnstructured(double     index,
                   double     diameter,
                   double     wavelength,
                   double     nMedium,
                   ndarray    Phi,
                   ndarray    Theta,
                   double     Polarization,
                   double     E0,
                   double     R)
{

  py::buffer_info PhiBuffer   = Phi.request();
  py::buffer_info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0];

  py::array_t<complex128> result0 = py::array_t<complex128>(Phi.size());
  py::array_t<complex128> result1 = py::array_t<complex128>(Phi.size());

  py::buffer_info buf0 = result0.request();
  py::buffer_info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
         k = 2. * PI / wavelength,
         temp0 ;

  complex128 *ptr0 = (complex128 *) buf0.ptr,
             *ptr1 = (complex128 *) buf1.ptr,
             temp2,
             propagator = E0 / (k * R) * exp(-j*k*R),
             *s1s2 = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  _S1S2(index, diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  for (auto k=0; k < PhiLength; k++ )
  {
    temp0 = ThetaPtr[k] ;
    ptr0[k] = j* propagator * s1s2[k] * (complex128) abs(cos(temp0 + Polarization));
    ptr1[k] =- propagator * s1s2[k + PhiLength] * (complex128) abs(sin(temp0 + Polarization));

  }

  free(s1s2);

  return std::make_pair(result0.attr("transpose")(), result1.attr("transpose")());
}










std::pair<Cndarray, Cndarray>
FieldsStructuredUnpolarized(double     index,
                             double     diameter,
                             double     wavelength,
                             double     nMedium,
                             ndarray    Phi,
                             ndarray    Theta,
                             double     E0,
                             double     R)
{

  py::buffer_info PhiBuffer   = Phi.request();
  py::buffer_info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0],
      w            = 0;

  py::array_t<complex128> result0 = py::array_t<complex128>(PhiLength * ThetaLength);
  py::array_t<complex128> result1 = py::array_t<complex128>(PhiLength * ThetaLength);

  py::buffer_info buf0 = result0.request();
  py::buffer_info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          k = 2. * PI / wavelength,
          temp0, temp1 = 1./sqrt(2.) ;

  complex128 *ptr0 = (complex128 *) buf0.ptr,
             *ptr1 = (complex128 *) buf1.ptr,
              propagator = E0 / (k * R) * exp(-j*k*R),
             *s1s2 = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  _S1S2(index, diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          temp0 = ThetaPtr[t] ;
          ptr0[w] = j* propagator * s1s2[p] * temp1;
          ptr1[w] =- propagator * s1s2[p + PhiLength] * temp1;
          w++;
       }
  }

  free(s1s2);
  result0.resize({PhiLength,ThetaLength});
  result1.resize({PhiLength,ThetaLength});

  return std::make_pair(result0.attr("transpose")(), result1.attr("transpose")());
}





std::pair<Cndarray, Cndarray>
FieldsUnstructuredUnpolarized(double     index,
                               double     diameter,
                               double     wavelength,
                               double     nMedium,
                               ndarray    Phi,
                               ndarray    Theta,
                               double     E0,
                               double     R)
{

  py::buffer_info PhiBuffer   = Phi.request();
  py::buffer_info ThetaBuffer = Theta.request();

  int PhiLength    = PhiBuffer.shape[0],
      ThetaLength  = ThetaBuffer.shape[0];

  py::array_t<complex128> result0 = py::array_t<complex128>(Phi.size());
  py::array_t<complex128> result1 = py::array_t<complex128>(Phi.size());

  py::buffer_info buf0 = result0.request();
  py::buffer_info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          k = 2. * PI / wavelength,
          temp0, temp1 = 1./sqrt(2.) ;

  complex128 *ptr0 = (complex128 *) buf0.ptr,
             *ptr1 = (complex128 *) buf1.ptr,
              propagator = E0 / (k * R) * exp(-j*k*R),
             *s1s2 = (complex128*) calloc(2 * PhiLength , sizeof(complex128));

  _S1S2(index, diameter, wavelength, nMedium, PhiPtr, PhiLength, s1s2);

  for (auto k=0; k < PhiLength; k++ )
  {
    temp0 = ThetaPtr[k] ;
    ptr0[k] = j* propagator * s1s2[k] * temp1;
    ptr1[k] =- propagator * s1s2[k + PhiLength] * temp1;
  }

  free(s1s2);

  return std::make_pair(result0.attr("transpose")(), result1.attr("transpose")());
}




}




// -
