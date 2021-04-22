#include <iostream>

class BASE{

    public:
        double E0, k, Polarization;

        void  ComputePrefactor(double* prefactor, uint MaxOrder);

        virtual void ComputeAnBn(complex128* an, complex128* bn, uint MaxOrder){};

        virtual double& Getk(){return this->k;};

        virtual double& GetE0(){return this->E0;};

        virtual double& GetPolarization(){return this->Polarization;};

        virtual double& GetSizeParam(){};

        std::tuple<double, double, double, double, double, double, double> GetEfficiencies();

        std::tuple<Cndarray,Cndarray> S1S2(const ndarray Phi),
                                      sS1S2(  ndarray& Phi, ndarray& Theta),
                                      uS1S2(  ndarray& Phi, ndarray& Theta),
                                      sFields(ndarray& Phi, ndarray& Theta, double R),
                                      uFields(ndarray& Phi, ndarray& Theta, double R);

        inline void MiePiTau(double mu, uint MaxOrder, complex128 *pin, complex128 *taun);

        BASE(){}

        virtual ~BASE(){ }

};

inline void
BASE::MiePiTau(double        mu,
               uint          MaxOrder,
               complex128   *pin,
               complex128   *taun)

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
BASE::ComputePrefactor(double* prefactor, uint MaxOrder)
{
   for (uint m = 0; m < MaxOrder ; m++)
   {
      prefactor[m] = (double) ( 2 * (m+1) + 1 ) / ( (m+1) * ( (m+1) + 1 ) );
   }
}



std::tuple<Cndarray,Cndarray>
BASE::S1S2(const ndarray Phi)
{
  uint MaxOrder           = GetMaxOrder(this->GetSizeParam());

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
BASE::sFields(ndarray& Phi, ndarray& Theta, double R)
{
  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0];

 Cndarray     ETheta     = Cndarray(PhiLength*ThetaLength),
              EPhi       = Cndarray(PhiLength*ThetaLength);

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * CosTerm      = (double*) calloc(ThetaLength, sizeof(double)),
             * SinTerm      = (double*) calloc(ThetaLength, sizeof(double));

  Cndarray     S1, S2;

  complex128   propagator = this->GetE0() / (this->Getk() * R) * exp(-JJ*this->Getk()*R);

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, this->GetPolarization());


  std::tie(S1, S2) = S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

   Structured(ThetaLength, PhiLength, S2Ptr, CosTerm, - propagator, EPhiPtr);

   Structured(ThetaLength, PhiLength, S1Ptr, SinTerm, JJ * propagator, EThetaPtr);

   EPhi.resize({PhiLength,ThetaLength});
   ETheta.resize({PhiLength,ThetaLength});

   EPhi   = EPhi.attr("transpose")();
   ETheta = ETheta.attr("transpose")();

   return std::make_tuple(EPhi, ETheta)  ;
}



std::tuple<Cndarray,Cndarray>
BASE::uFields(ndarray& Phi, ndarray& Theta, double R)
{
  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0];

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * CosTerm      = (double*) calloc(ThetaLength, sizeof(double)),
             * SinTerm      = (double*) calloc(ThetaLength, sizeof(double));

  Cndarray     ETheta     = Cndarray(PhiLength),
               EPhi       = Cndarray(PhiLength),
               S1,
               S2;

  complex128   propagator = this->GetE0() / (this->Getk() * R) * exp(-JJ*this->Getk()*R);

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, this->GetPolarization());

  std::tie(S1, S2) = this->S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

  Unstructured(ThetaLength, PhiLength, S2Ptr, CosTerm, - propagator, EPhiPtr);

  Unstructured(ThetaLength, PhiLength, S1Ptr, SinTerm, JJ * propagator, EThetaPtr);

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(EPhi, ETheta)  ;
}




std::tuple<Cndarray,Cndarray>
BASE::sS1S2(ndarray& Phi, ndarray& Theta)
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

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, this->GetPolarization());

  std::tie(S1, S2) = this->S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

  Structured(ThetaLength, PhiLength, S2Ptr, CosTerm, this->GetE0(), EPhiPtr);

  Structured(ThetaLength, PhiLength, S1Ptr, SinTerm, this->GetE0(), EThetaPtr);

  EPhi.resize({PhiLength,ThetaLength});
  ETheta.resize({PhiLength,ThetaLength});

  EPhi   = EPhi.attr("transpose")();
  ETheta = ETheta.attr("transpose")();

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(EPhi, ETheta)  ;
}



std::tuple<Cndarray,Cndarray>
BASE::uS1S2(ndarray& Phi, ndarray& Theta)
{

  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0];

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * CosTerm      = (double*) calloc(ThetaLength, sizeof(double)),
             * SinTerm      = (double*) calloc(ThetaLength, sizeof(double));

  Cndarray     ETheta     = Cndarray(PhiLength),
               EPhi       = Cndarray(PhiLength),
               S1,
               S2;

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, this->GetPolarization());

  std::tie(S1, S2) = this->S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

  Unstructured(ThetaLength, PhiLength, S2Ptr, CosTerm, this->GetE0(), EPhiPtr);

  Unstructured(ThetaLength, PhiLength, S1Ptr, SinTerm, this->GetE0(), EThetaPtr);

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(EPhi, ETheta)  ;

}




std::tuple<double, double, double, double, double, double, double>
BASE::GetEfficiencies()
{
    uint MaxOrder   = GetMaxOrder(this->GetSizeParam());

    complex128 * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
               * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128));

    this->ComputeAnBn(an, bn, MaxOrder);

    double Qsca   = GetQsca(an, bn, MaxOrder, this->GetSizeParam());
    double Qext   = GetQext(an, bn, MaxOrder, this->GetSizeParam());
    double g      = Getg(an, bn, MaxOrder, this->GetSizeParam(), Qsca);
    double Qabs   = Qext - Qsca;
    double Qback  = GetQback(an, bn, MaxOrder, this->GetSizeParam());
    double Qpr    = Qext - g * Qsca;
    double Qratio = Qback / Qsca;

    free(an);
    free(bn);
    return std::make_tuple(Qsca, Qext, Qabs, Qback, Qratio, g, Qpr);
}

// -
