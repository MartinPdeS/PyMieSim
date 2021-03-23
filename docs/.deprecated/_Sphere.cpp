#include <iostream>

#include <math.h>


class SPHERE{

private:
      bool       Polarized;

      double     Diameter,
                 Index,
                 nMedium,
                 SizeParam,
                 Polarization,
                 Wavelength,
                 k,
                 E0,
                 Qsca,
                 Qext,
                 Qabs;

      int        MaxOrder;


      std::tuple<Cndarray, Cndarray>  PrivateGetS1S2(ndarray Phi);

      complex128 *_an,
                 *_bn,
                 *pin,
                 *taun,
                  propagator,
                 *An_data,
                 *Bn_data,
                 *Cn_data,
                 *Dn_data,
                 *S1_data,
                 *S2_data;

        void      GetQext(),
                  GetQsca(),
                  CoefficientAnBn(),
                  LowFrequencyMie_ab(),
                  HighFrequencyMie_ab(),
                  MiePiTau(double mu),

                  LoopStructuredPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr, complex128 *EPhiPtr, complex128 *EThetaPtr, Cndarray  S1, Cndarray  S2),

                  LoopUnstructuredPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr, complex128 *EPhiPtr, complex128 *EThetaPtr),

                  LoopStructuredUPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr, complex128 *EPhiPtr, complex128 *EThetaPtr),

                  LoopUnstructuredUPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr, complex128 *EPhiPtr, complex128 *EThetaPtr);

    public:
        Cndarray&                          PublicAn(int MaxOrder);
        Cndarray&                          PublicBn(int MaxOrder);
        Cndarray&                          PublicCn(int MaxOrder);
        Cndarray&                          PublicDn(int MaxOrder);
        std::tuple<Cndarray, Cndarray>     FieldsStructured(ndarray Phi, ndarray Theta, double R);
        std::tuple<Cndarray&, Cndarray&>   FieldsUnstructured(ndarray Phi, ndarray Theta, double R);
        std::tuple<Cndarray&, Cndarray&>   PublicGetS1S2 (ndarray Phi);
        std::tuple<double, double, double> GetEfficiencies();


  SPHERE(double Index,
         double Diameter,
         double Wavelength,
         double nMedium,
         double Polarization,
         double E0)
    {
        Diameter      = Diameter;
        Index         = Index;
        nMedium       = nMedium;
        Wavelength    = Wavelength;
        E0            = E0;
        k             = 2 * PI / Wavelength;
        Polarization  = Polarization;
        SizeParam     = GetSizeParameter(Diameter, Wavelength, nMedium);
        MaxOrder      = GetMaxOrder(SizeParam);

        _an           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        _bn           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        pin           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        taun          = (complex128*) malloc(sizeof(complex128)*MaxOrder);

        Cndarray An         = Cndarray(MaxOrder);
        Cndarray Bn         = Cndarray(MaxOrder);
        Cndarray Cn         = Cndarray(MaxOrder);
        Cndarray Dn         = Cndarray(MaxOrder);

        info AnBuffer       = An.request();
        info BnBuffer       = Bn.request();
        info CnBuffer       = Cn.request();
        info DnBuffer       = Dn.request();

        An_data       = (complex128 *) AnBuffer.ptr;
        Bn_data       = (complex128 *) BnBuffer.ptr;
        Cn_data       = (complex128 *) CnBuffer.ptr;
        Dn_data       = (complex128 *) DnBuffer.ptr;

        if (Polarization == -1.){ Polarized = false  ; }
        else                    { Polarized = true   ; }
    }
};



void
SPHERE::CoefficientAnBn()
{
  if (SizeParam < 0.5){LowFrequencyMie_ab(); }

  else{HighFrequencyMie_ab();}
}


void
SPHERE::LowFrequencyMie_ab()
{
  double LL, m2, x3, x4, x5, x6;

  m2          = Index * Index;
  LL          = (m2 - 1) / (m2 + 2);
  x3          = SizeParam * SizeParam * SizeParam;
  x4          = x3 * SizeParam;
  x5          = x4 * SizeParam;
  x6          = x5 * SizeParam;

  _an[0] = (-2.*J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  _an[1] = (-1.*J * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  _bn[0] = (-1.*J * x5 / 45.) * (m2 - 1.);
  _bn[1] = 0. + 0.*J;

  MaxOrder = 2;
}


void
SPHERE::HighFrequencyMie_ab()
{

  const double mx = Index * SizeParam,
               temp  = sqrt(0.5 * PI * SizeParam);

  const int nmx = (int) ( std::max( MaxOrder, (int) abs(mx) ) + 16 );

  iVec gsx, gs1x, px, chx, p1x, ch1x, D, da, db;

  Vec Dn = Vec(nmx);

  p1x.push_back( sin(SizeParam) );
  ch1x.push_back( cos(SizeParam) );

  for (double i = nmx - 1; i > 1; i--){  Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );}

  for (auto i = 0; i < MaxOrder; i++)
    {
        px.push_back(  temp * Jn( (double)(i+1)+0.5, SizeParam ) );
        chx.push_back(-temp * Yn( (double)(i+1)+0.5, SizeParam ) );

        p1x.push_back(px[i]);
        ch1x.push_back(chx[i]);

        gsx.push_back( px[i] - 1.*J * chx[i] );
        gs1x.push_back( p1x[i] - 1.*J * ch1x[i] );

        da.push_back( Dn[i+1] / Index + (double)(i+1) / SizeParam );
        db.push_back( Index * Dn[i+1] + (double)(i+1) / SizeParam );

        _an[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;
        _bn[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;
    }
}


void
SPHERE::GetQsca()
{
    complex128 temp     = 0.;

    for(auto it = 0; it < MaxOrder; ++it)
    {
         temp += (2.* (double)(it+1) + 1.) * (   std::real( _an[it] ) * std::real( _an[it] )
                                               + std::imag( _an[it] ) * std::imag( _an[it] )
                                               + std::real( _bn[it] ) * std::real( _bn[it] )
                                               + std::imag( _bn[it] ) * std::imag( _bn[it] ) );
    }
    Qsca = 2. / (SizeParam * SizeParam)  * std::real(temp);
}


void
SPHERE::GetQext()
{
    complex128 temp     = 0.;
    for(auto it = 0; it < MaxOrder; ++it)
    {
      temp += ( 2.*(double)(it+1) + 1.) * ( std::real( _an[it] + _an[it] ) );
    }

    Qext = 2. / (SizeParam * SizeParam) * std::real(temp);
}


std::tuple<double, double, double>
SPHERE::GetEfficiencies()
{
    CoefficientAnBn();
    GetQsca();
    GetQext();
    Qabs = Qext - Qsca;

    return std::make_tuple(Qsca, Qext, Qabs);
}


void
SPHERE::MiePiTau(double mu)

{
  pin[0] = 1.;
  pin[1] = 3. * mu;

  taun[0] = mu;
  taun[1] = 3.0 * cos(2. * acos(mu) );

  double n = 0;
  for (auto i = 2; i < MaxOrder; i++)
      {
       n = (double)i;

       pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

       taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
     }
}


std::tuple<Cndarray, Cndarray>
SPHERE::PrivateGetS1S2(ndarray Phi)
{
    info        PhiBuffer = Phi.request();

    double     *PhiPtr    = (double *) PhiBuffer.ptr,
                prefactor = 0.;

    int         PhiLength = PhiBuffer.size;

    Cndarray S1           = Cndarray(PhiLength);
    Cndarray S2           = Cndarray(PhiLength);

    info S1Buffer         = S1.request(),
         S2Buffer         = S2.request();

    S1_data       = (complex128 *) S1Buffer.ptr;
    S2_data       = (complex128 *) S2Buffer.ptr;

    complex128  temp0     = 0., temp1   = 0.;

    CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < MaxOrder ; k++){
            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( _an[k] * pin[k] +  _bn[k] * taun[k] );
            temp1    += prefactor * ( _an[k] * taun[k] + _bn[k] * pin[k]  );
          }

        S1_data[i] = temp0;
        S2_data[i] = temp1;

        temp0 = 0.; temp1=0.;
    }
    free(pin); free(taun);
    return std::tie<Cndarray, Cndarray>(S1, S2);
}



std::tuple<Cndarray&, Cndarray&>
SPHERE::PublicGetS1S2(ndarray Phi)
{
    info        PhiBuffer = Phi.request();

    double     *PhiPtr    = (double *) PhiBuffer.ptr,
                prefactor = 0.;

    int         PhiLength = PhiBuffer.size;

    complex128  temp0        = 0.,
                temp1        = 0.;

    Cndarray S1            = Cndarray(PhiLength);
    Cndarray S2            = Cndarray(PhiLength);

    info S1Buffer       = S1.request(),
         S2Buffer       = S2.request();

    S1_data       = (complex128 *) S1Buffer.ptr;
    S2_data       = (complex128 *) S2Buffer.ptr;

    CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < MaxOrder ; k++){

            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( _an[k] * pin[k] +  _bn[k] * taun[k] );
            temp1    += prefactor * ( _an[k] * taun[k] + _bn[k] * pin[k]  );
          }

        S1_data[i] = temp0;
        S2_data[i] = temp1;


        temp0 = 0.; temp1=0.;

    }
    //free(pin); free(taun);
    return std::tuple<Cndarray&, Cndarray&>(S1, S2);
}


Cndarray&
SPHERE::PublicAn(int MaxOrder)
{
  const double mx          = Index * SizeParam,
               temp        = sqrt(0.5 * PI * SizeParam);

  const int    nmx         = (int) ( std::max( MaxOrder, (int) abs(mx) ) + 16 );

  iVec gsx, gs1x, px, chx, p1x, ch1x, D, da, db;

  Vec Dn = Vec(nmx);

  p1x.push_back( sin(SizeParam) );
  ch1x.push_back( cos(SizeParam) );

  for (double i = nmx - 1; i > 1; i--){  Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );}

  for (auto i = 0; i < MaxOrder; i++)
    {
        px.push_back(  temp * Jn( (double)(i+1)+0.5, SizeParam ) );
        chx.push_back(-temp * Yn( (double)(i+1)+0.5, SizeParam ) );

        p1x.push_back(px[i]);
        ch1x.push_back(chx[i]);

        gsx.push_back(   px[i] - 1.*J * chx[i]  );
        gs1x.push_back( p1x[i] - 1.*J * ch1x[i] );

        da.push_back( Dn[i+1] / Index + (double)(i+1) / SizeParam );

        An_data[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;

    }
  //return An;
}



Cndarray&
SPHERE::PublicBn(int MaxOrder)
{
  const double mx          = Index * SizeParam,
               temp        = sqrt(0.5 * PI * SizeParam);

  const int    nmx         = (int) ( std::max( MaxOrder, (int) abs(mx) ) + 16 );

  iVec gsx, gs1x, px, chx, p1x, ch1x, D, da, db;

  Vec Dn = Vec(nmx);

  p1x.push_back( sin(SizeParam) );
  ch1x.push_back( cos(SizeParam) );

  for (double i = nmx - 1; i > 1; i--){  Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );}

  for (auto i = 0; i < MaxOrder; i++)
    {
        px.push_back(  temp * Jn( (double)(i+1)+0.5, SizeParam ) );
        chx.push_back(-temp * Yn( (double)(i+1)+0.5, SizeParam ) );

        p1x.push_back(px[i]);
        ch1x.push_back(chx[i]);

        gsx.push_back(   px[i] - 1.*J * chx[i]  );
        gs1x.push_back( p1x[i] - 1.*J * ch1x[i] );

        db.push_back( Index * Dn[i+1] + (double)(i+1) / SizeParam );

        Bn_data[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;

    }
  //return Bn;
}



Cndarray&
SPHERE::PublicCn(int MaxOrder)
{
  const double mx          = Index * SizeParam,
               x           = SizeParam,
               temp        = sqrt(0.5 * PI * SizeParam),
               MuSp        = 1.,
               Mu          = 1.,
               M           = Index/nMedium;


  complex128 numerator,
             denominator;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator        = M * MuSp * ( Xi(order, x) * Psi_p(order, x) - Xi_p(order, x) * Psi(order, x) );
    denominator      = MuSp * Xi(order, x) * Psi_p(order, mx) - Mu * M * Xi_p(order, x) * Psi(order, mx);
    Cn_data[order-1] = numerator / denominator;
  }
  //return Cn;
}



Cndarray&
SPHERE::PublicDn(int MaxOrder)
{
  const double mx          = Index * SizeParam,
               x           = SizeParam,
               temp        = sqrt(0.5 * PI * SizeParam),
               MuSp        = 1.,
               Mu          = 1.,
               M           = Index/nMedium;

  complex128 numerator,
             denominator;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator        = Mu * M*M * ( Xi(order, x) * Psi_p(order, x) - Xi_p(order, x) * Psi(order, x) );
    denominator      = Mu * M * Xi(order, x) * Psi_p(order, mx) - MuSp * Xi_p(order, x) * Psi(order, mx);
    Dn_data[order-1] = numerator / denominator;
  }
  //return Dn;
}






// -
