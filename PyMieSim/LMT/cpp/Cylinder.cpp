#include <iostream>
#pragma GCC visibility push(hidden)





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
                 Mu,
                 MuScat,
                 GetQsca(),
                 GetQabs(),
                 GetQext();

    void        ComputePrefactor(double* prefactor, uint MaxOrder),
                ComputeAnBn(complex128* anPtr, complex128* bnPtr, uint MaxOrder),
                LowFreqAnBn(complex128* anPtr, complex128* bnPtr),
                HighFreqAnBn(complex128* anPtr, complex128* bnPtr, uint MaxOrder),
                Structured(uint ThetaLength, uint PhiLength, complex128 *array0, double *array1, complex128  scalar, complex128 *output),
                Unstructured(uint ThetaLength, uint PhiLength, complex128 *array0, double *array1, complex128  scalar, complex128 *output);


    inline void MiePiTau(double mu, uint MaxOrder, complex128 *pin, complex128 *taun);



    public:
        std::tuple<double, double, double> GetEfficiencies();
        Cndarray                           An(uint MaxOrder);
        Cndarray                           Bn(uint MaxOrder);
        Cndarray                           Cn(uint MaxOrder);
        Cndarray                           Dn(uint MaxOrder);
        std::tuple<Cndarray,Cndarray>      S1S2(ndarray Phi);
        std::tuple<Cndarray,Cndarray>      SFields(ndarray Phi, ndarray Theta, double R);
        std::tuple<Cndarray,Cndarray>      UFields(ndarray Phi, ndarray Theta, double R);



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
        this->Mu            = 1.0;
        this->MuScat        = 1.0;
        }

        ~CYLINDER(){  }

};



void
CYLINDER::ComputePrefactor(double* prefactor, uint MaxOrder)
{
   for (uint m = 0; m < MaxOrder ; m++)
   {
      prefactor[m] = (double) ( 2 * (m+1) + 1 ) / ( (m+1) * ( (m+1) + 1 ) );
   }
}


inline void
CYLINDER::MiePiTau(double      mu,
                 uint        MaxOrder,
                 complex128 *pin,
                 complex128 *taun)

{
  pin[0] = 1.;
  pin[1] = 3. * mu;

  taun[0] = mu;
  taun[1] = 3.0 * cos(2. * acos(mu) );

  double n = 0;
  for (uint i = 2; i < MaxOrder; i++)
      {
       n = (double)i;

       pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

       taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
     }
}



void
CYLINDER::ComputeAnBn(complex128* anPtr, complex128* bnPtr, uint MaxOrder)
{
  return HighFreqAnBn(anPtr, bnPtr, MaxOrder) ;
}

void
CYLINDER::HighFreqAnBn(complex128* anPtr, complex128* bnPtr, uint MaxOrder)
{

  double x = SizeParam;

  complex128 numerator, denominator;

  for (int order = 1; order < MaxOrder+1; order++)
  {
      numerator   = Index * Jn(order, Index*x) * Jn_p(order, nMedium*x) - nMedium * Jn_p(order, Index*x) * Jn(order, nMedium*x);
      denominator = Index * Jn(order, Index*x) * Hn_p(order, nMedium*x) - nMedium * Jn_p(order, Index*x) * Hn(order, nMedium*x);
      anPtr[order-1] = numerator/denominator;

      numerator   = nMedium * Jn(order, Index*x) * Jn_p(order, nMedium*x) - Index*Jn_p(order, Index*x) * Jn(order, nMedium*x);
      denominator = nMedium * Jn(order, Index*x) * Hn_p(order, nMedium*x) - Index*Jn_p(order, Index*x) * Hn(order, nMedium*x);
      bnPtr[order-1] = numerator/denominator;
  }
}


void
CYLINDER::LowFreqAnBn(complex128* anPtr, complex128* bnPtr)
{

  double LL, m2, x3, x4, x5, x6;

  m2          = Index * Index;
  LL          = (m2 - 1) / (m2 + 2);
  x3          = SizeParam * SizeParam * SizeParam;
  x4          = x3 * SizeParam;
  x5          = x4 * SizeParam;
  x6          = x5 * SizeParam;

  anPtr[0] = (-2.*J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  anPtr[1] = (-1.*J * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  bnPtr[0] = (-1.*J * x5 / 45.) * (m2 - 1.);
  bnPtr[1] = 0. + 0.*J;

}



std::tuple<Cndarray,Cndarray>
CYLINDER::S1S2(const ndarray Phi)
{

  uint MaxOrder           = GetMaxOrder(SizeParam);

  uint PhiLength          = Phi.request().shape[0];

  double     * PhiPtr     = (double*) Phi.request().ptr,
             * prefactor  = (double*) calloc(MaxOrder, sizeof(double));

  Cndarray     s1         = Cndarray(PhiLength),
               s2         = Cndarray(PhiLength);

  complex128 * s1Ptr      = (complex128 *) s1.request().ptr,
             * s2Ptr      = (complex128 *) s2.request().ptr,
             * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
             * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
             * pin        = (complex128*) calloc(MaxOrder, sizeof(complex128)),
             * taun       = (complex128*) calloc(MaxOrder, sizeof(complex128));

  this->ComputeAnBn(an, bn, MaxOrder);

  this->ComputePrefactor(prefactor, MaxOrder);

  for (uint i = 0; i < PhiLength; i++){

      MiePiTau( cos( PhiPtr[i]-PI/2 ), MaxOrder, pin, taun );
      s1Ptr[i] = 0.;
      s2Ptr[i] = 0.;

      for (uint m = 0; m < MaxOrder ; m++){
          s1Ptr[i]    += prefactor[m] * ( an[m] * pin[m] +  bn[m] * taun[m] );
          s2Ptr[i]    += prefactor[m] * ( an[m] * taun[m] + bn[m] * pin[m]  );
        }
  }

  free(pin);
  free(taun);
  free(an);
  free(bn);
  free(prefactor);
  s1Ptr = NULL;
  s2Ptr = NULL;

  return std::make_tuple(s1, s2)  ;
}


std::tuple<Cndarray,Cndarray>
CYLINDER::SFields(ndarray Phi, ndarray Theta, double R)
{

  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0];

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * CosTerm      = (double*) calloc(ThetaLength, sizeof(double)),
             * SinTerm      = (double*) calloc(ThetaLength, sizeof(double));

  Cndarray     ETheta     = Cndarray(PhiLength*ThetaLength),
               EPhi       = Cndarray(PhiLength*ThetaLength),
               S1,
               S2;

  complex128   propagator = E0 / (k * R) * exp(-J*k*R);

  for (uint t = 0; t < ThetaLength; t++)
  {
    CosTerm[t] = abs(cos(Polarization + ThetaPtr[t])) ;
    SinTerm[t] = abs(sin(Polarization + ThetaPtr[t])) ;
  }

  std::tie(S1, S2) = this->S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

  Structured(ThetaLength, PhiLength, S2Ptr, SinTerm, - propagator, EThetaPtr);

  Structured(ThetaLength, PhiLength, S1Ptr, CosTerm, J * propagator, EPhiPtr);

  EPhi.resize({PhiLength,ThetaLength});
  ETheta.resize({PhiLength,ThetaLength});

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(ETheta, EPhi)  ;

}

void
CYLINDER::Structured(uint         ThetaLength,
                   uint         PhiLength,
                   complex128 *array0,
                   double     *array1,
                   complex128  scalar,
                   complex128 *output)
{
  for (uint p=0; p < PhiLength; p++ )
  {
    for (uint t=0; t < ThetaLength; t++ )
    {
      *output   = scalar * array0[p] * array1[t];
       output++;
    }
  }
}


void
CYLINDER::Unstructured(uint        ThetaLength,
                     uint        PhiLength,
                     complex128 *array0,
                     double     *array1,
                     complex128  scalar,
                     complex128 *output)
{
  for (uint p=0; p < PhiLength; p++ )
  {
    *output   = scalar * array0[p] * array1[p];
     output++;
  }
}


std::tuple<Cndarray,Cndarray>
CYLINDER::UFields(ndarray Phi, ndarray Theta, double R)
{

  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0];

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * CosTerm      = (double*) calloc(ThetaLength, sizeof(double)),
             * SinTerm      = (double*) calloc(ThetaLength, sizeof(double));

  Cndarray     ETheta     = Cndarray(PhiLength*ThetaLength),
               EPhi       = Cndarray(PhiLength*ThetaLength),
               S1,
               S2;

  complex128   propagator = E0 / (k * R) * exp(-J*k*R);

  for (uint t = 0; t < ThetaLength; t++)
  {
    CosTerm[t] = abs(cos(Polarization + ThetaPtr[t])) ;
    SinTerm[t] = abs(sin(Polarization + ThetaPtr[t])) ;
  }

  std::tie(S1, S2) = this->S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

  Unstructured(ThetaLength, PhiLength, S2Ptr, SinTerm, - propagator, EThetaPtr);

  Unstructured(ThetaLength, PhiLength, S1Ptr, CosTerm, J * propagator, EPhiPtr);

  EPhi.resize({PhiLength,ThetaLength});
  ETheta.resize({PhiLength,ThetaLength});

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(ETheta, EPhi)  ;

}


Cndarray
CYLINDER::Bn(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return bn;
}


Cndarray
CYLINDER::An(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return an;
}



double
CYLINDER::GetQsca()
{
    uint MaxOrder   = GetMaxOrder(SizeParam);

    complex128 * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
               * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128));

    this->ComputeAnBn(an, bn, MaxOrder);

    complex128 temp = 0.;

    for(uint it = 0; it < MaxOrder; ++it)
    {
         temp += (2.* (double)(it+1) + 1.) * (   std::real( an[it] ) * std::real( an[it] )
                                               + std::imag( an[it] ) * std::imag( an[it] )
                                               + std::real( bn[it] ) * std::real( bn[it] )
                                               + std::imag( bn[it] ) * std::imag( bn[it] ) );
    }

    free(an);
    free(bn);

    return 2. / (SizeParam * SizeParam)  * std::real(temp);
}


double
CYLINDER::GetQext()
{
    uint MaxOrder   = GetMaxOrder(SizeParam);

    complex128 * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
               * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128));

    this->ComputeAnBn(an, bn, MaxOrder);

    complex128 temp = 0.;

    for(uint it = 0; it < MaxOrder; ++it){ temp += ( 2.*(double)(it+1) + 1.) * ( std::real( an[it] + an[it] ) ); }

    free(an);
    free(bn);

    return 2. / (SizeParam * SizeParam) * std::real(temp);
}


std::tuple<double, double, double>
CYLINDER::GetEfficiencies()
{
    double Qsca = GetQsca();
    double Qext = GetQext();
    double Qabs = Qext - Qsca;

    return std::make_tuple(Qsca, Qext, Qabs);
}
