#include <iostream>




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

      Cndarray    An,
                  Bn,
                  Cn,
                  Dn,
                  S1,
                  S2,
                  EPhi,
                  ETheta;

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
                 *S2_data,
                 *EPhi_data,
                 *ETheta_data;

        void      GetQext(),
                  GetQsca(),
                  CoefficientAnBn(),
                  LowFrequencyMie_ab(),
                  HighFrequencyMie_ab(),
                  MiePiTau(double mu),
                  PrivateGetS1S2(ndarray Phi),

                  LoopStructuredPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr),

                  LoopUnstructuredPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr),

                  LoopStructuredUPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr),

                  LoopUnstructuredUPol(int PhiLength, int ThetaLength, double *PhiPtr, double *ThetaPtr);

    public:
        Cndarray&                          PublicAn(int MaxOrder);
        Cndarray&                          PublicBn(int MaxOrder);
        Cndarray&                          PublicCn(int MaxOrder);
        Cndarray&                          PublicDn(int MaxOrder);
        std::tuple<Cndarray&, Cndarray&>   FieldsStructured(ndarray Phi, ndarray Theta, double R);
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
        this->Diameter      = Diameter;
        this->Index         = Index;
        this->nMedium       = nMedium;
        this->Wavelength    = Wavelength;
        this->E0            = E0;
        this->k             = 2 * PI / Wavelength;
        this->Polarization  = Polarization;
        this->SizeParam     = GetSizeParameter(Diameter, Wavelength, nMedium);
        this->MaxOrder      = GetMaxOrder(this->SizeParam);

        this->_an           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->_bn           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->pin           = (complex128*) malloc(sizeof(complex128)*MaxOrder);
        this->taun          = (complex128*) malloc(sizeof(complex128)*MaxOrder);

        this->An            = Cndarray(MaxOrder);
        this->Bn            = Cndarray(MaxOrder);
        this->Cn            = Cndarray(MaxOrder);
        this->Dn            = Cndarray(MaxOrder);

        info AnBuffer       = this->An.request(),
             BnBuffer       = this->Bn.request(),
             CnBuffer       = this->Cn.request(),
             DnBuffer       = this->Dn.request();

        this->An_data       = (complex128 *) AnBuffer.ptr;
        this->Bn_data       = (complex128 *) BnBuffer.ptr;
        this->Cn_data       = (complex128 *) CnBuffer.ptr;
        this->Dn_data       = (complex128 *) DnBuffer.ptr;

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

  this->_an[0] = (-2.*J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  this->_an[1] = (-1.*J * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  this->_bn[0] = (-1.*J * x5 / 45.) * (m2 - 1.);
  this->_bn[1] = 0. + 0.*J;

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

        this->_an[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;
        this->_bn[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;
    }
}


void
SPHERE::GetQsca()
{
    complex128 temp     = 0.;

    for(auto it = 0; it < this->MaxOrder; ++it)
    {
         temp += (2.* (double)(it+1) + 1.) * (   std::real( this->_an[it] ) * std::real( this->_an[it] )
                                               + std::imag( this->_an[it] ) * std::imag( this->_an[it] )
                                               + std::real( this->_bn[it] ) * std::real( this->_bn[it] )
                                               + std::imag( this->_bn[it] ) * std::imag( this->_bn[it] ) );
    }
    this->Qsca = 2. / (this->SizeParam * this->SizeParam)  * std::real(temp);
}


void
SPHERE::GetQext()
{
    complex128 temp     = 0.;
    for(auto it = 0; it < this->MaxOrder; ++it)
    {
      temp += ( 2.*(double)(it+1) + 1.) * ( std::real( this->_an[it] + this->_an[it] ) );
    }

    this->Qext = 2. / (this->SizeParam * this->SizeParam) * std::real(temp);
}


std::tuple<double, double, double>
SPHERE::GetEfficiencies()
{
    this->CoefficientAnBn();
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
SPHERE::PrivateGetS1S2(ndarray Phi)
{
    info        PhiBuffer = Phi.request();

    double     *PhiPtr    = (double *) PhiBuffer.ptr,
                prefactor = 0.;

    int         PhiLength = PhiBuffer.size;

    this->S1            = Cndarray(PhiLength);
    this->S2            = Cndarray(PhiLength);

    info S1Buffer       = this->S1.request(),
         S2Buffer       = this->S2.request();

    this->S1_data       = (complex128 *) S1Buffer.ptr;
    this->S2_data       = (complex128 *) S2Buffer.ptr;

    complex128  temp0     = 0., temp1   = 0.;

    this->CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        this->MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < this->MaxOrder ; k++){
            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( this->_an[k] * this->pin[k] +  this->_bn[k] * this->taun[k] );
            temp1    += prefactor * ( this->_an[k] * this->taun[k] + this->_bn[k] * this->pin[k]  );
          }

        this->S1_data[i] = temp0;
        this->S2_data[i] = temp1;

        temp0 = 0.; temp1=0.;
    }
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

    this->S1            = Cndarray(PhiLength);
    this->S2            = Cndarray(PhiLength);

    info S1Buffer       = this->S1.request(),
         S2Buffer       = this->S2.request();

    this->S1_data       = (complex128 *) S1Buffer.ptr;
    this->S2_data       = (complex128 *) S2Buffer.ptr;

    this->CoefficientAnBn();

    for (auto i = 0; i < PhiLength; i++){

        this->MiePiTau( cos( PhiPtr[i]-PI/2 ) );

        for (auto k = 0; k < this->MaxOrder ; k++){

            prefactor = (double) ( 2 * (k+1) + 1 ) / ( (k+1) * ( (k+1) + 1 ) );
            temp0    += prefactor * ( this->_an[k] * this->pin[k] +  this->_bn[k] * this->taun[k] );
            temp1    += prefactor * ( this->_an[k] * this->taun[k] + this->_bn[k] * this->pin[k]  );
          }

        this->S1_data[i] = temp0;
        this->S2_data[i] = temp1;


        temp0 = 0.; temp1=0.;
    }
    return std::tuple<Cndarray&, Cndarray&>(this->S1, this->S2);
}


Cndarray&
SPHERE::PublicAn(int MaxOrder)
{
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

        this->An_data[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;

    }
  return this->An;
}



Cndarray&
SPHERE::PublicBn(int MaxOrder)
{
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

        this->Bn_data[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;

    }
  return this->Bn;
}



Cndarray&
SPHERE::PublicCn(int MaxOrder)
{
  const double mx          = this->Index * this->SizeParam,
               x           = this->SizeParam,
               temp        = sqrt(0.5 * PI * this->SizeParam),
               MuSp        = 1.,
               Mu          = 1.,
               M           = this->Index/this->nMedium;


  complex128 numerator,
             denominator;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator        = M * MuSp * ( Xi(order, x) * Psi_p(order, x) - Xi_p(order, x) * Psi(order, x) );
    denominator      = MuSp * Xi(order, x) * Psi_p(order, mx) - Mu * M * Xi_p(order, x) * Psi(order, mx);
    this->Cn_data[order-1] = numerator / denominator;
  }
  return this->Cn;
}



Cndarray&
SPHERE::PublicDn(int MaxOrder)
{
  const double mx          = this->Index * this->SizeParam,
               x           = this->SizeParam,
               temp        = sqrt(0.5 * PI * this->SizeParam),
               MuSp        = 1.,
               Mu          = 1.,
               M           = this->Index/this->nMedium;

  complex128 numerator,
             denominator;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator        = Mu * M*M * ( Xi(order, x) * Psi_p(order, x) - Xi_p(order, x) * Psi(order, x) );
    denominator      = Mu * M * Xi(order, x) * Psi_p(order, mx) - MuSp * Xi_p(order, x) * Psi(order, mx);
    this->Dn_data[order-1] = numerator / denominator;
  }
  return this->Dn;
}




void
SPHERE::LoopStructuredPol(int         PhiLength,
                          int         ThetaLength,
                          double     *PhiPtr,
                          double     *ThetaPtr)
{
  int w = 0; double temp0;
  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          temp0          = ThetaPtr[t] ;
          this->EPhi_data[w]   = J * this->propagator * S1_data[p] * (complex128) abs(cos(temp0 + this->Polarization));
          this->ETheta_data[w] = - this->propagator * S2_data[p] * (complex128) abs(sin(temp0 + this->Polarization));
          w++;
       }
  }
}


void
SPHERE::LoopStructuredUPol(int         PhiLength,
                           int         ThetaLength,
                           double     *PhiPtr,
                           double     *ThetaPtr)
{
  int w = 0;
  double temp0 = 1./sqrt(2.);
  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          this->EPhi_data[w] = J * this->propagator * this->S1_data[p] * temp0;
          this->ETheta_data[w] = - this->propagator * this->S2_data[p] * temp0;
          w++;
       }
  }
}



void
SPHERE::LoopUnstructuredPol(int         PhiLength,
                            int         ThetaLength,
                            double     *PhiPtr,
                            double     *ThetaPtr)
{
  double temp0;
    for (auto p=0; p < PhiLength; p++ )
    {
        temp0 = ThetaPtr[p] ;
        this->EPhi_data[p]   = J * this->propagator * this->S1_data[p] * (complex128) abs(cos(temp0 + this->Polarization));
        this->ETheta_data[p] =   - this->propagator * this->S2_data[p] * (complex128) abs(sin(temp0 + this->Polarization));
    }
}


void
SPHERE::LoopUnstructuredUPol(int        PhiLength,
                            int         ThetaLength,
                            double     *PhiPtr,
                            double     *ThetaPtr)
  {
  double temp0 = 1./sqrt(2.);
    for (auto p=0; p < PhiLength; p++ )
    {
        this->EPhi_data[p]   = J * this->propagator * this->S1_data[p] * temp0 ;
        this->ETheta_data[p] =   - this->propagator * this->S2_data[p] * temp0 ;
    }
}


std::tuple<Cndarray&, Cndarray&>
SPHERE::FieldsStructured(ndarray Phi, ndarray Theta, double R)
{
  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0];

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr;

  this->EPhi          = Cndarray(PhiLength*ThetaLength);
  this->ETheta        = Cndarray(PhiLength*ThetaLength);

  info EPhiBuffer          = this->EPhi.request(),
       EThetaBuffer        = this->ETheta.request();

  this->EPhi_data     = (complex128 *) EPhiBuffer.ptr;
  this->ETheta_data   = (complex128 *) EThetaBuffer.ptr;

  this->propagator   = this->E0 / (this->k * R) * exp(-J*this->k*R);

  this->PrivateGetS1S2(Phi);

  if (this->Polarized){ LoopStructuredPol (PhiLength, ThetaLength, PhiPtr, ThetaPtr); }
  else                { LoopStructuredUPol(PhiLength, ThetaLength, PhiPtr, ThetaPtr); }

  this->EPhi.resize({PhiLength,ThetaLength});//.attr("transpose")();
  this->ETheta.resize({PhiLength,ThetaLength});//.attr("transpose")();



  return std::tie<Cndarray&, Cndarray&>(this->EPhi, this->ETheta);

}


std::tuple<Cndarray&, Cndarray&>
SPHERE::FieldsUnstructured(ndarray Phi, ndarray Theta, double R)
{
  info        PhiBuffer    = Phi.request(),
              ThetaBuffer  = Theta.request();

  int         PhiLength    = PhiBuffer.shape[0],
              ThetaLength  = ThetaBuffer.shape[0];

  double     *PhiPtr       = (double *) PhiBuffer.ptr,
             *ThetaPtr     = (double *) ThetaBuffer.ptr;

  this->EPhi               = Cndarray(PhiLength);
  this->ETheta             = Cndarray(PhiLength);

  info EPhiBuffer          = this->EPhi.request(),
       EThetaBuffer        = this->ETheta.request();

  this->EPhi_data          = (complex128 *) EPhiBuffer.ptr;
  this->ETheta_data        = (complex128 *) EThetaBuffer.ptr;

  this->propagator   = this->E0 / (this->k * R) * exp(-J*this->k*R);

  this->PrivateGetS1S2(Phi);

  if (this->Polarized){ LoopUnstructuredPol (PhiLength, ThetaLength, PhiPtr, ThetaPtr); }
  else                { LoopUnstructuredUPol(PhiLength, ThetaLength, PhiPtr, ThetaPtr); }

  return std::tie<Cndarray&, Cndarray&>(this->EPhi, this->ETheta);

}








// -
