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


  CYLINDER(double Index,
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

    this->_an[order-1] = numerator/denominator;

    numerator   = m * Jn(order, mt*x) * Jn_p(order, m*x) - mt*Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = m * Jn(order, mt*x) * Hn_p(order, m*x) - mt*Jn_p(order, mt*x) * Hn(order, m*x);
    this->_bn[order-1] = numerator/denominator;
  }
}


void
CYLINDER::GetQsca()
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
CYLINDER::GetQext()
{
    complex128 temp     = 0.;
    for(auto it = 0; it < this->MaxOrder; ++it)
    {
      temp += ( 2.*(double)(it+1) + 1.) * ( std::real( this->_an[it] + this->_an[it] ) );
    }

    this->Qext = 2. / (this->SizeParam * this->SizeParam) * std::real(temp);
}


std::tuple<double, double, double>
CYLINDER::GetEfficiencies()
{
    this->CoefficientAnBn();
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
CYLINDER::PrivateGetS1S2(ndarray Phi)
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

    complex128  temp0     = 0.,
                temp1   = 0.;

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



Cndarray&
CYLINDER::PublicAn(int MaxOrder)
{
  complex128 numerator,
             denominator;

  double x  = this->SizeParam,
         m  = this->nMedium,
         mt = this->Index;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator   = mt * Jn(order, mt*x) * Jn_p(order, m*x) - m * Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = mt * Jn(order, mt*x) * Hn_p(order, m*x) - m * Jn_p(order, mt*x) * Hn(order, m*x);
    An_data[order-1] = numerator/denominator;
  }
  return this->An;
}



Cndarray&
CYLINDER::PublicBn(int MaxOrder)
{
  complex128 numerator,
             denominator;

  double x  = this->SizeParam,
         m  = this->nMedium,
         mt = this->Index;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator   = m * Jn(order, mt*x) * Jn_p(order, m*x) - mt*Jn_p(order, mt*x) * Jn(order, m*x);
    denominator = m * Jn(order, mt*x) * Hn_p(order, m*x) - mt*Jn_p(order, mt*x) * Hn(order, m*x);
    this->Bn_data[order-1] = numerator/denominator;
  }
  return this->Bn;
}


std::tuple<Cndarray&, Cndarray&>
CYLINDER::PublicGetS1S2(ndarray Phi)
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


void
CYLINDER::LoopStructuredPol(int         PhiLength,
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
CYLINDER::LoopStructuredUPol(int         PhiLength,
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
CYLINDER::LoopUnstructuredPol(int         PhiLength,
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
CYLINDER::LoopUnstructuredUPol(int        PhiLength,
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
CYLINDER::FieldsStructured(ndarray Phi, ndarray Theta, double R)
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
CYLINDER::FieldsUnstructured(ndarray Phi, ndarray Theta, double R)
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
