#include "Math.cpp"
#include "utils.cpp"
#include <iostream>





class CYLINDER{

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

  complex128 *an,
             *bn,
             *pin,
             *taun,
              propagator;

public:
  void        CoefficientAnBn(),
              LowFrequencyMie_ab(),
              HighFrequencyMie_ab(),
              GetQsca(),
              GetQext(),
              MiePiTau(double mu),
              PrivateGetS1S2(ndarray Phi, complex128* S1_data, complex128* S2_data),

              LoopStructuredPol(int PhiLength,
                                int ThetaLength,
                                double *PhiPtr,
                                double *ThetaPtr,
                                complex128 *EPhi_data,
                                complex128 *ETheta_data,
                                complex128 *S1_data,
                                complex128 *S2_data),

              LoopUnstructuredPol(int PhiLength,
                                  int ThetaLength,
                                  double *PhiPtr,
                                  double *ThetaPtr,
                                  complex128 *EPhi_data,
                                  complex128 *ETheta_data,
                                  complex128 *S1_data,
                                  complex128 *S2_data),

              LoopStructuredUPol(int PhiLength,
                                 int ThetaLength,
                                 double *PhiPtr,
                                 double *ThetaPtr,
                                 complex128 *EPhi_data,
                                 complex128 *ETheta_data,
                                 complex128 *S1_data,
                                 complex128 *S2_data),

              LoopUnstructuredUPol(int PhiLength,
                                   int ThetaLength,
                                   double *PhiPtr,
                                   double *ThetaPtr,
                                   complex128 *EPhi_data,
                                   complex128 *ETheta_data,
                                   complex128 *S1_data,
                                   complex128 *S2_data);

  Cndarray                           PublicAn(int MaxOrder);
  Cndarray                           PublicBn(int MaxOrder);
  std::tuple<Cndarray, Cndarray>     FieldsStructured(ndarray Phi, ndarray Theta, double R);
  std::tuple<Cndarray, Cndarray>     FieldsUnstructured(ndarray Phi, ndarray Theta, double R);
  std::tuple<Cndarray, Cndarray>     PublicGetS1S2 (ndarray Phi);
  std::tuple<double, double, double> GetEfficiencies();


  CYLINDER(double Index,
           double Diameter,
           double Wavelength,
           double nMedium,
           double Polarization,
           double E0)
    {
        this->Diameter     = Diameter;
        this->Index        = Index;
        this->nMedium      = nMedium;
        this->Wavelength   = Wavelength;
        this->E0           = E0;
        this->k            = 2 * PI / Wavelength;
        this->Polarization = Polarization;
        this->SizeParam    = GetSizeParameter(Diameter, Wavelength, nMedium);
        this->MaxOrder     = GetMaxOrder(this->SizeParam);
        this->an           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->bn           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->pin          = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->taun         = (complex128*) malloc(sizeof(complex128)*MaxOrder);

        if (Polarization == -1.){ this->Polarized = false  ; }
        else                    { this->Polarized = true   ; }
    }
};



void
CYLINDER::CoefficientAnBn() // ref: https://doi.org/10.1364/AO.44.002338
{
  complex128 numerator, denominator;
  double x  = this->SizeParam,
         m  = this->nMedium,
         mt = this->Index;

  for (auto order = 1; order < this->MaxOrder+1; order++)
  {

    numerator   = mt * Jn(order, mt*x) * Jn_p(order, m*x) - m * Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = mt * Jn(order, mt*x) * Hn_p(order, m*x) - m * Jn_p(order, mt*x) * Hn(order, m*x);

    this->an[order-1] = numerator/denominator;

    numerator   = m * Jn(order, mt*x) * Jn_p(order, m*x) - mt*Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = m * Jn(order, mt*x) * Hn_p(order, m*x) - mt*Jn_p(order, mt*x) * Hn(order, m*x);
    this->bn[order-1] = numerator/denominator;
  }
}


void
CYLINDER::GetQsca()
{
    complex128 temp     = 0.;

    for(auto it = 0; it < this->MaxOrder; ++it)
    {
         temp += (2.* (double)(it+1) + 1.) * (   std::real( this->an[it] ) * std::real( this->an[it] )
                                               + std::imag( this->an[it] ) * std::imag( this->an[it] )
                                               + std::real( this->bn[it] ) * std::real( this->bn[it] )
                                               + std::imag( this->bn[it] ) * std::imag( this->bn[it] ) );
    }
    this->Qsca = 2. / (this->SizeParam * this->SizeParam)  * std::real(temp);
}


void
CYLINDER::GetQext()
{
    complex128 temp     = 0.;
    for(auto it = 0; it < this->MaxOrder; ++it)
    {
      temp += ( 2.*(double)(it+1) + 1.) * ( std::real( this->an[it] + this->an[it] ) );
    }

    this->Qext = 2. / (this->SizeParam * this->SizeParam) * std::real(temp);
}


std::tuple<double, double, double>
CYLINDER::GetEfficiencies()
{
    this->GetQsca();
    this->GetQext();
    this->Qabs = this->Qext - Qsca;

    return std::make_tuple(this->Qsca, this->Qext, this->Qabs);
}


void
CYLINDER::MiePiTau(double mu)

{
  this->pin[0] = 1.;
  this->pin[1] = 3. * mu;

  this->taun[0] = mu;
  this->taun[1] = 3.0 * cos(2. * acos(mu) );

  double n = 0;
  for (auto i = 2; i < this->MaxOrder; i++)
      {
       n = (double)i;

       this->pin[i] = ( (2. * n + 1.) * mu * this->pin[i-1] - (n + 1.) * this->pin[i-2] ) / n;

       this->taun[i] = (n + 1.) * mu * this->pin[i] - (n + 2.) * this->pin[i-1];
     }
}


void
CYLINDER::PrivateGetS1S2(ndarray Phi, complex128* S1_data, complex128* S2_data)
{
    info        PhiBuffer = Phi.request();

    double     *PhiPtr    = (double *) PhiBuffer.ptr,
                prefactor = 0.;

    int         PhiLength = PhiBuffer.size;

    complex128  temp0     = 0., temp1   = 0.;

    this->CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        this->MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < this->MaxOrder ; k++){

            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( this->an[k] * this->pin[k] +  this->bn[k] * this->taun[k] );
            temp1    += prefactor * ( this->an[k] * this->taun[k] + this->bn[k] * this->pin[k]  );
          }

        S1_data[i] = temp0;
        S2_data[i] = temp1;


        temp0 = 0.; temp1=0.;
    }
}



Cndarray
CYLINDER::PublicAn(int MaxOrder)
{
  complex128 numerator,
             denominator;

  double x  = this->SizeParam,
         m  = this->nMedium,
         mt = this->Index;

  Cndarray    An           = Cndarray(MaxOrder);

  info        AnBuf        = An.request();

  complex128 *An_data      = (complex128 *) AnBuf.ptr;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator   = mt * Jn(order, mt*x) * Jn_p(order, m*x) - m * Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = mt * Jn(order, mt*x) * Hn_p(order, m*x) - m * Jn_p(order, mt*x) * Hn(order, m*x);
    An_data[order-1] = numerator/denominator;
  }
  return An;
}



Cndarray
CYLINDER::PublicBn(int MaxOrder)
{
  complex128 numerator,
             denominator;

  double x  = this->SizeParam,
         m  = this->nMedium,
         mt = this->Index;

  Cndarray    Bn           = Cndarray(MaxOrder);

  info        BnBuf        = Bn.request();

  complex128 *Bn_data      = (complex128 *) BnBuf.ptr;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator   = m * Jn(order, mt*x) * Jn_p(order, m*x) - mt*Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = m * Jn(order, mt*x) * Hn_p(order, m*x) - mt*Jn_p(order, mt*x) * Hn(order, m*x);
    Bn_data[order-1] = numerator/denominator;
  }
  return Bn;
}


std::tuple<Cndarray, Cndarray>
CYLINDER::PublicGetS1S2(ndarray Phi)
{
    info        PhiBuffer    = Phi.request();

    double     *PhiPtr       = (double *) PhiBuffer.ptr,
                prefactor    = 0.;

    int         PhiLength    = PhiBuffer.size;

    Cndarray    S1           = Cndarray(PhiLength),
                S2           = Cndarray(PhiLength);

    info        S1Buf        = S1.request(),
                S2Buf        = S2.request();

    complex128 *S1_data      = (complex128 *) S1Buf.ptr,
               *S2_data      = (complex128 *) S2Buf.ptr,
                temp0        = 0.,
                temp1        = 0.;

    this->CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        this->MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < this->MaxOrder ; k++){

            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( this->an[k] * this->pin[k] +  this->bn[k] * this->taun[k] );
            temp1    += prefactor * ( this->an[k] * this->taun[k] + this->bn[k] * this->pin[k]  );
          }

        S1_data[i] = temp0;
        S2_data[i] = temp1;


        temp0 = 0.; temp1=0.;
    }
    return std::tie<Cndarray, Cndarray>(S1, S2);
}


void
CYLINDER::LoopStructuredPol(int         PhiLength,
                          int         ThetaLength,
                          double     *PhiPtr,
                          double     *ThetaPtr,
                          complex128 *EPhi_data,
                          complex128 *ETheta_data,
                          complex128 *S1_data,
                          complex128 *S2_data)
{
  int w = 0; double temp0;
  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          temp0 = ThetaPtr[t] ;
          EPhi_data[w] = J * this->propagator * S1_data[p] * (complex128) abs(cos(temp0 + this->Polarization));
          ETheta_data[w] = - this->propagator * S2_data[p] * (complex128) abs(sin(temp0 + this->Polarization));
          w++;
       }
  }
}


void
CYLINDER::LoopStructuredUPol(int         PhiLength,
                           int         ThetaLength,
                           double     *PhiPtr,
                           double     *ThetaPtr,
                           complex128 *EPhi_data,
                           complex128 *ETheta_data,
                           complex128 *S1_data,
                           complex128 *S2_data)
{
  int w = 0;
  double temp0 = 1./sqrt(2.);
  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          EPhi_data[w] = J * this->propagator * S1_data[p] * temp0;
          ETheta_data[w] = - this->propagator * S2_data[p] * temp0;
          w++;
       }
  }
}



void
CYLINDER::LoopUnstructuredPol(int         PhiLength,
                            int         ThetaLength,
                            double     *PhiPtr,
                            double     *ThetaPtr,
                            complex128 *EPhi_data,
                            complex128 *ETheta_data,
                            complex128 *S1_data,
                            complex128 *S2_data)
{
  double temp0;
    for (auto p=0; p < PhiLength; p++ )
    {
        temp0 = ThetaPtr[p] ;
        EPhi_data[p]   = J * this->propagator * S1_data[p] * (complex128) abs(cos(temp0 + this->Polarization));
        ETheta_data[p] =   - this->propagator * S2_data[p] * (complex128) abs(sin(temp0 + this->Polarization));
    }
}


void
CYLINDER::LoopUnstructuredUPol(int         PhiLength,
                            int         ThetaLength,
                            double     *PhiPtr,
                            double     *ThetaPtr,
                            complex128 *EPhi_data,
                            complex128 *ETheta_data,
                            complex128 *S1_data,
                            complex128 *S2_data)
  {
  double temp0 = 1./sqrt(2.);
    for (auto p=0; p < PhiLength; p++ )
    {
        EPhi_data[p]   = J * this->propagator * S1_data[p] * temp0 ;
        ETheta_data[p] =   - this->propagator * S2_data[p] * temp0 ;
    }
}


std::tuple<Cndarray, Cndarray>
CYLINDER::FieldsStructured(ndarray Phi, ndarray Theta, double R)
{
  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0];

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr;

  Cndarray    EPhi         = Cndarray(PhiLength * ThetaLength),
              ETheta       = Cndarray(PhiLength * ThetaLength),
              S1           = Cndarray(PhiLength),
              S2           = Cndarray(PhiLength);

  info        EPhiBuf         = EPhi.request(),
              EThetaBuf       = ETheta.request(),
              S1Buf           = S1.request(),
              S2Buf           = S2.request();

  complex128 *EPhi_data    = (complex128 *) EPhiBuf.ptr,
             *ETheta_data  = (complex128 *) EThetaBuf.ptr,
             *S1_data      = (complex128 *) S1Buf.ptr,
             *S2_data      = (complex128 *) S2Buf.ptr;

  this->propagator   = this->E0 / (this->k * R) * exp(-J*this->k*R);

  this->PrivateGetS1S2(Phi, S1_data, S2_data);

  if (this->Polarized){ LoopStructuredPol (PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }
  else                { LoopStructuredUPol(PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }

  EPhi.resize({PhiLength,ThetaLength});
  ETheta.resize({PhiLength,ThetaLength});

  return std::tuple<Cndarray, Cndarray>(EPhi.attr("transpose")(), ETheta.attr("transpose")());

}


std::tuple<Cndarray, Cndarray>
CYLINDER::FieldsUnstructured(ndarray Phi, ndarray Theta, double R)
{
  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0];

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr;

  Cndarray    EPhi         = Cndarray(PhiLength * ThetaLength),
              ETheta       = Cndarray(PhiLength * ThetaLength),
              S1           = Cndarray(PhiLength),
              S2           = Cndarray(PhiLength);

  info        EPhiBuf         = EPhi.request(),
              EThetaBuf       = ETheta.request(),
              S1Buf           = S1.request(),
              S2Buf           = S2.request();

  complex128 *EPhi_data    = (complex128 *) EPhiBuf.ptr,
             *ETheta_data  = (complex128 *) EThetaBuf.ptr,
             *S1_data      = (complex128 *) S1Buf.ptr,
             *S2_data      = (complex128 *) S2Buf.ptr;

  this->propagator   = this->E0 / (this->k * R) * exp(-J*this->k*R);

  this->PrivateGetS1S2(Phi, S1_data, S2_data);

  if (this->Polarized){ LoopUnstructuredPol (PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }
  else                { LoopUnstructuredUPol(PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }

  return std::tie<Cndarray, Cndarray>(EPhi, ETheta);

}






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

  complex128 *an,
             *bn,
             *pin,
             *taun,
              propagator;

public:
  void        CoefficientAnBn(),
              LowFrequencyMie_ab(),
              HighFrequencyMie_ab(),
              GetQsca(),
              GetQext(),
              MiePiTau(double mu),
              PrivateGetS1S2(ndarray Phi, complex128* S1_data, complex128* S2_data),

              LoopStructuredPol(int PhiLength,
                                int ThetaLength,
                                double *PhiPtr,
                                double *ThetaPtr,
                                complex128 *EPhi_data,
                                complex128 *ETheta_data,
                                complex128 *S1_data,
                                complex128 *S2_data),

              LoopUnstructuredPol(int PhiLength,
                                  int ThetaLength,
                                  double *PhiPtr,
                                  double *ThetaPtr,
                                  complex128 *EPhi_data,
                                  complex128 *ETheta_data,
                                  complex128 *S1_data,
                                  complex128 *S2_data),

              LoopStructuredUPol(int PhiLength,
                                 int ThetaLength,
                                 double *PhiPtr,
                                 double *ThetaPtr,
                                 complex128 *EPhi_data,
                                 complex128 *ETheta_data,
                                 complex128 *S1_data,
                                 complex128 *S2_data),

              LoopUnstructuredUPol(int PhiLength,
                                   int ThetaLength,
                                   double *PhiPtr,
                                   double *ThetaPtr,
                                   complex128 *EPhi_data,
                                   complex128 *ETheta_data,
                                   complex128 *S1_data,
                                   complex128 *S2_data);

  Cndarray                           PublicAn(int MaxOrder);
  Cndarray                           PublicBn(int MaxOrder);
  std::tuple<Cndarray, Cndarray>     FieldsStructured(ndarray Phi, ndarray Theta, double R);
  std::tuple<Cndarray, Cndarray>     FieldsUnstructured(ndarray Phi, ndarray Theta, double R);
  std::tuple<Cndarray, Cndarray>     PublicGetS1S2 (ndarray Phi);
  std::tuple<double, double, double> GetEfficiencies();


  SPHERE(double Index,
         double Diameter,
         double Wavelength,
         double nMedium,
         double Polarization,
         double E0)
    {
        this->Diameter     = Diameter;
        this->Index        = Index;
        this->nMedium      = nMedium;
        this->Wavelength   = Wavelength;
        this->E0           = E0;
        this->k            = 2 * PI / Wavelength;
        this->Polarization = Polarization;
        this->SizeParam    = GetSizeParameter(Diameter, Wavelength, nMedium);
        this->MaxOrder     = GetMaxOrder(this->SizeParam);
        this->an           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->bn           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->pin          = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->taun         = (complex128*) malloc(sizeof(complex128)*MaxOrder);

        if (Polarization == -1.){ this->Polarized = false  ; }
        else                    { this->Polarized = true   ; }
    }
};

void
SPHERE::CoefficientAnBn()
{
  if (this->SizeParam < 0.5){this->LowFrequencyMie_ab(); }

  else{this->HighFrequencyMie_ab();}
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

  this->an[0] = (-2.*J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  this->an[1] = (-1.*J * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  this->bn[0] = (-1.*J * x5 / 45.) * (m2 - 1.);
  this->bn[1] = 0. + 0.*J;

  this->MaxOrder = 2;
}


void
SPHERE::HighFrequencyMie_ab()
{

  const double mx = Index * this->SizeParam,
               temp  = sqrt(0.5 * PI * this->SizeParam);

  const int nmx = (int) ( std::max( this->MaxOrder, (int) abs(mx) ) + 16 );

  iVec gsx, gs1x, px, chx, p1x, ch1x, D, da, db;

  Vec Dn = Vec(nmx);

  p1x.push_back( sin(this->SizeParam) );
  ch1x.push_back( cos(this->SizeParam) );

  for (double i = nmx - 1; i > 1; i--){  Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );}

  for (auto i = 0; i < this->MaxOrder; i++)
    {
        px.push_back(  temp * Jn( (double)(i+1)+0.5, this->SizeParam ) );
        chx.push_back(-temp * Yn( (double)(i+1)+0.5, this->SizeParam ) );

        p1x.push_back(px[i]);
        ch1x.push_back(chx[i]);

        gsx.push_back( px[i] - 1.*J * chx[i] );
        gs1x.push_back( p1x[i] - 1.*J * ch1x[i] );

        da.push_back( Dn[i+1] / Index + (double)(i+1) / this->SizeParam );
        db.push_back( Index * Dn[i+1] + (double)(i+1) / this->SizeParam );

        this->an[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;
        this->bn[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;
    }
}


void
SPHERE::GetQsca()
{
    complex128 temp     = 0.;

    for(auto it = 0; it < this->MaxOrder; ++it)
    {
         temp += (2.* (double)(it+1) + 1.) * (   std::real( this->an[it] ) * std::real( this->an[it] )
                                               + std::imag( this->an[it] ) * std::imag( this->an[it] )
                                               + std::real( this->bn[it] ) * std::real( this->bn[it] )
                                               + std::imag( this->bn[it] ) * std::imag( this->bn[it] ) );
    }
    this->Qsca = 2. / (this->SizeParam * this->SizeParam)  * std::real(temp);
}


void
SPHERE::GetQext()
{
    complex128 temp     = 0.;
    for(auto it = 0; it < this->MaxOrder; ++it)
    {
      temp += ( 2.*(double)(it+1) + 1.) * ( std::real( this->an[it] + this->an[it] ) );
    }

    this->Qext = 2. / (this->SizeParam * this->SizeParam) * std::real(temp);
}


std::tuple<double, double, double>
SPHERE::GetEfficiencies()
{
    this->GetQsca();
    this->GetQext();
    this->Qabs = this->Qext - Qsca;

    return std::make_tuple(this->Qsca, this->Qext, this->Qabs);
}


void
SPHERE::MiePiTau(double mu)

{
  this->pin[0] = 1.;
  this->pin[1] = 3. * mu;

  this->taun[0] = mu;
  this->taun[1] = 3.0 * cos(2. * acos(mu) );

  double n = 0;
  for (auto i = 2; i < this->MaxOrder; i++)
      {
       n = (double)i;

       this->pin[i] = ( (2. * n + 1.) * mu * this->pin[i-1] - (n + 1.) * this->pin[i-2] ) / n;

       this->taun[i] = (n + 1.) * mu * this->pin[i] - (n + 2.) * this->pin[i-1];
     }
}


void
SPHERE::PrivateGetS1S2(ndarray Phi, complex128* S1_data, complex128* S2_data)
{
    info        PhiBuffer = Phi.request();

    double     *PhiPtr    = (double *) PhiBuffer.ptr,
                prefactor = 0.;

    int         PhiLength = PhiBuffer.size;

    complex128  temp0     = 0., temp1   = 0.;

    this->CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        this->MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < this->MaxOrder ; k++){

            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( this->an[k] * this->pin[k] +  this->bn[k] * this->taun[k] );
            temp1    += prefactor * ( this->an[k] * this->taun[k] + this->bn[k] * this->pin[k]  );
          }

        S1_data[i] = temp0;
        S2_data[i] = temp1;


        temp0 = 0.; temp1=0.;
    }
}



std::tuple<Cndarray, Cndarray>
SPHERE::PublicGetS1S2(ndarray Phi)
{
    info        PhiBuffer = Phi.request();

    double     *PhiPtr    = (double *) PhiBuffer.ptr,
                prefactor = 0.;

    int         PhiLength = PhiBuffer.size;

    Cndarray    S1           = Cndarray(PhiLength),
                S2           = Cndarray(PhiLength);

    info        S1Buf           = S1.request(),
                S2Buf           = S2.request();

    complex128 *S1_data      = (complex128 *) S1Buf.ptr,
               *S2_data      = (complex128 *) S2Buf.ptr,
                temp0        = 0.,
                temp1        = 0.;

    this->CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        this->MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < this->MaxOrder ; k++){

            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( this->an[k] * this->pin[k] +  this->bn[k] * this->taun[k] );
            temp1    += prefactor * ( this->an[k] * this->taun[k] + this->bn[k] * this->pin[k]  );
          }

        S1_data[i] = temp0;
        S2_data[i] = temp1;


        temp0 = 0.; temp1=0.;
    }
    return std::tie<Cndarray, Cndarray>(S1, S2);
}


Cndarray
SPHERE::PublicAn(int MaxOrder)
{
  Cndarray     An          = Cndarray(MaxOrder);

  info         AnBuf       = An.request();

  complex128  *An_data     = (complex128 *) AnBuf.ptr;

  const double mx          = Index * this->SizeParam,
               temp        = sqrt(0.5 * PI * this->SizeParam);

  const int    nmx         = (int) ( std::max( this->MaxOrder, (int) abs(mx) ) + 16 );

  iVec gsx, gs1x, px, chx, p1x, ch1x, D, da, db;

  Vec Dn = Vec(nmx);

  p1x.push_back( sin(this->SizeParam) );
  ch1x.push_back( cos(this->SizeParam) );

  for (double i = nmx - 1; i > 1; i--){  Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );}

  for (auto i = 0; i < this->MaxOrder; i++)
    {
        px.push_back(  temp * Jn( (double)(i+1)+0.5, this->SizeParam ) );
        chx.push_back(-temp * Yn( (double)(i+1)+0.5, this->SizeParam ) );

        p1x.push_back(px[i]);
        ch1x.push_back(chx[i]);

        gsx.push_back(   px[i] - 1.*J * chx[i]  );
        gs1x.push_back( p1x[i] - 1.*J * ch1x[i] );

        da.push_back( Dn[i+1] / Index + (double)(i+1) / this->SizeParam );

        An_data[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;

    }
  return An;
}



Cndarray
SPHERE::PublicBn(int MaxOrder)
{
  Cndarray     Bn          = Cndarray(MaxOrder);

  info         BnBuf       = Bn.request();

  complex128  *Bn_data     = (complex128 *) BnBuf.ptr;

  const double mx          = Index * this->SizeParam,
               temp        = sqrt(0.5 * PI * this->SizeParam);

  const int    nmx         = (int) ( std::max( this->MaxOrder, (int) abs(mx) ) + 16 );

  iVec gsx, gs1x, px, chx, p1x, ch1x, D, da, db;

  Vec Dn = Vec(nmx);

  p1x.push_back( sin(this->SizeParam) );
  ch1x.push_back( cos(this->SizeParam) );

  for (double i = nmx - 1; i > 1; i--){  Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );}

  for (auto i = 0; i < this->MaxOrder; i++)
    {
        px.push_back(  temp * Jn( (double)(i+1)+0.5, this->SizeParam ) );
        chx.push_back(-temp * Yn( (double)(i+1)+0.5, this->SizeParam ) );

        p1x.push_back(px[i]);
        ch1x.push_back(chx[i]);

        gsx.push_back(   px[i] - 1.*J * chx[i]  );
        gs1x.push_back( p1x[i] - 1.*J * ch1x[i] );

        db.push_back( Index * Dn[i+1] + (double)(i+1) / this->SizeParam );

        Bn_data[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;

    }
  return Bn;
}



void
SPHERE::LoopStructuredPol(int         PhiLength,
                          int         ThetaLength,
                          double     *PhiPtr,
                          double     *ThetaPtr,
                          complex128 *EPhi_data,
                          complex128 *ETheta_data,
                          complex128 *S1_data,
                          complex128 *S2_data)
{
  int w = 0; double temp0;
  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          temp0 = ThetaPtr[t] ;
          EPhi_data[w] = J * this->propagator * S1_data[p] * (complex128) abs(cos(temp0 + this->Polarization));
          ETheta_data[w] = - this->propagator * S2_data[p] * (complex128) abs(sin(temp0 + this->Polarization));
          w++;
       }
  }
}


void
SPHERE::LoopStructuredUPol(int         PhiLength,
                           int         ThetaLength,
                           double     *PhiPtr,
                           double     *ThetaPtr,
                           complex128 *EPhi_data,
                           complex128 *ETheta_data,
                           complex128 *S1_data,
                           complex128 *S2_data)
{
  int w = 0;
  double temp0 = 1./sqrt(2.);
  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          EPhi_data[w] = J * this->propagator * S1_data[p] * temp0;
          ETheta_data[w] = - this->propagator * S2_data[p] * temp0;
          w++;
       }
  }
}



void
SPHERE::LoopUnstructuredPol(int         PhiLength,
                            int         ThetaLength,
                            double     *PhiPtr,
                            double     *ThetaPtr,
                            complex128 *EPhi_data,
                            complex128 *ETheta_data,
                            complex128 *S1_data,
                            complex128 *S2_data)
{
  double temp0;
    for (auto p=0; p < PhiLength; p++ )
    {
        temp0 = ThetaPtr[p] ;
        EPhi_data[p]   = J * this->propagator * S1_data[p] * (complex128) abs(cos(temp0 + this->Polarization));
        ETheta_data[p] =   - this->propagator * S2_data[p] * (complex128) abs(sin(temp0 + this->Polarization));
    }
}


void
SPHERE::LoopUnstructuredUPol(int         PhiLength,
                            int         ThetaLength,
                            double     *PhiPtr,
                            double     *ThetaPtr,
                            complex128 *EPhi_data,
                            complex128 *ETheta_data,
                            complex128 *S1_data,
                            complex128 *S2_data)
  {
  double temp0 = 1./sqrt(2.);
    for (auto p=0; p < PhiLength; p++ )
    {
        EPhi_data[p]   = J * this->propagator * S1_data[p] * temp0 ;
        ETheta_data[p] =   - this->propagator * S2_data[p] * temp0 ;
    }
}


std::tuple<Cndarray, Cndarray>
SPHERE::FieldsStructured(ndarray Phi, ndarray Theta, double R)
{
  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0];

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr;

  Cndarray    EPhi         = Cndarray(PhiLength * ThetaLength),
              ETheta       = Cndarray(PhiLength * ThetaLength),
              S1           = Cndarray(PhiLength),
              S2           = Cndarray(PhiLength);

  info        EPhiBuf      = EPhi.request(),
              EThetaBuf    = ETheta.request(),
              S1Buf        = S1.request(),
              S2Buf        = S2.request();

  complex128 *EPhi_data    = (complex128 *) EPhiBuf.ptr,
             *ETheta_data  = (complex128 *) EThetaBuf.ptr,
             *S1_data      = (complex128 *) S1Buf.ptr,
             *S2_data      = (complex128 *) S2Buf.ptr;

  this->propagator   = this->E0 / (this->k * R) * exp(-J*this->k*R);

  this->PrivateGetS1S2(Phi, S1_data, S2_data);

  if (this->Polarized){ LoopStructuredPol (PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }
  else                { LoopStructuredUPol(PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }

  EPhi.resize({PhiLength,ThetaLength});
  ETheta.resize({PhiLength,ThetaLength});

  return std::tuple<Cndarray, Cndarray>(EPhi.attr("transpose")(), ETheta.attr("transpose")());

}


std::tuple<Cndarray, Cndarray>
SPHERE::FieldsUnstructured(ndarray Phi, ndarray Theta, double R)
{
  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0];

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr;

  Cndarray    EPhi         = Cndarray(PhiLength),
              ETheta       = Cndarray(PhiLength),
              S1           = Cndarray(PhiLength),
              S2           = Cndarray(PhiLength);

  info        EPhiBuf      = EPhi.request(),
              EThetaBuf    = ETheta.request(),
              S1Buf        = S1.request(),
              S2Buf        = S2.request();

  complex128 *EPhi_data    = (complex128 *) EPhiBuf.ptr,
             *ETheta_data  = (complex128 *) EThetaBuf.ptr,
             *S1_data      = (complex128 *) S1Buf.ptr,
             *S2_data      = (complex128 *) S2Buf.ptr;

  this->propagator   = this->E0 / (this->k * R) * exp(-J*this->k*R);

  this->PrivateGetS1S2(Phi, S1_data, S2_data);

  if (this->Polarized){ LoopUnstructuredPol (PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }
  else                { LoopUnstructuredUPol(PhiLength, ThetaLength, PhiPtr, ThetaPtr, EPhi_data, ETheta_data, S1_data, S2_data); }

  return std::tie<Cndarray, Cndarray>(EPhi, ETheta);

}





PYBIND11_MODULE(Scatterer, module) {
    module.doc() = "LGeneralized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";


      py::class_<SPHERE>(module, "SPHERE")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("Index"),
           py::arg("Diameter"),
           py::arg("Wavelength"),
           py::arg("nMedium")      = 1.,
           py::arg("Polarization") = 0.,
           py::arg("E0")           = 1. )

      .def("S1S2",
           &SPHERE::PublicGetS1S2,
           py::arg("Phi")  )

      .def("SFields",
           &SPHERE::FieldsStructured,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R")  )

      .def("UFields",
           &SPHERE::FieldsUnstructured,
           py::arg("Phi"),
           py::arg("Theta"),
           py::arg("R")  )

      .def("an", &SPHERE::PublicAn, py::arg("MaxOrder")  = 5)

      .def("bn", &SPHERE::PublicBn, py::arg("MaxOrder")  = 5)

      .def_property("Efficiencies", &SPHERE::GetEfficiencies, &SPHERE::GetEfficiencies);


     py::class_<CYLINDER>(module, "CYLINDER")
     .def(py::init<double, double, double, double, double, double>(),
          py::arg("Index"),
          py::arg("Diameter"),
          py::arg("Wavelength"),
          py::arg("nMedium")      = 1.,
          py::arg("Polarization") = 0.,
          py::arg("E0")           = 1. )

     .def("S1S2",
          &CYLINDER::PublicGetS1S2,
          py::arg("Phi")  )

     .def("SFields",
          &CYLINDER::FieldsStructured,
          py::arg("Phi"),
          py::arg("Theta"),
          py::arg("R")  )

     .def("UFields",
          &CYLINDER::FieldsUnstructured,
          py::arg("Phi"),
          py::arg("Theta"),
          py::arg("R")  )

      .def("an", &CYLINDER::PublicAn, py::arg("MaxOrder")  = 5)

      .def("bn", &CYLINDER::PublicBn, py::arg("MaxOrder")  = 5)

      .def_property("Efficiencies", &CYLINDER::GetEfficiencies, &CYLINDER::GetEfficiencies);



}







// -
