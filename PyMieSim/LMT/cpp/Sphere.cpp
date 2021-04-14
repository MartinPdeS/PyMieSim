#include <iostream>
#include <math.h>




class SPHERE: public BASE{

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
             MuScat;

    void     ComputeAnBn(complex128* an, complex128* bn, uint MaxOrder),
             LowFreqAnBn(complex128* an, complex128* bn),
             HighFreqAnBn(complex128* an, complex128* bn, uint MaxOrder);


    public:
      std::tuple<double, double, double, double, double, double, double> GetEfficiencies();

      Cndarray                           An(uint MaxOrder),
                                         Bn(uint MaxOrder),
                                         Cn(uint MaxOrder),
                                         Dn(uint MaxOrder);

      std::tuple<Cndarray,Cndarray>      S1S2(ndarray Phi),
                                         sS1S2(ndarray& Phi, ndarray& Theta),
                                         uS1S2(ndarray& Phi, ndarray& Theta),
                                         sFields(ndarray& Phi, ndarray& Theta, double R),
                                         uFields(ndarray& Phi, ndarray& Theta, double R);



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
          this->Mu            = 1.0;
          this->MuScat        = 1.0;
        }

        ~SPHERE(){  }

};


void
SPHERE::ComputeAnBn(complex128* an, complex128* bn, uint MaxOrder)
{
  if (SizeParam < 0.5){LowFreqAnBn(an, bn) ; }
  else                {HighFreqAnBn(an, bn, MaxOrder) ; }
}


void
SPHERE::HighFreqAnBn(complex128* an, complex128* bn, uint MaxOrder)
{

  const double mx = Index * SizeParam, temp = sqrt(0.5 * PI * SizeParam);

  uint nmx = std::max( MaxOrder, (uint) mx ) + 16;

  iVec gsx, gs1x, px, chx, p1x, ch1x, D, da, db;

  Vec Dn = Vec(nmx);

  p1x.push_back( sin(SizeParam) );
  ch1x.push_back( cos(SizeParam) );

  for (double i = nmx - 1; i > 1; i--){  Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );}

  for (uint i = 0; i < MaxOrder; i++)
    {
        px.push_back(  temp * Jn( (double)(i+1)+0.5, SizeParam ) );
        chx.push_back(-temp * Yn( (double)(i+1)+0.5, SizeParam ) );

        p1x.push_back(px[i]);
        ch1x.push_back(chx[i]);

        gsx.push_back( px[i] - 1.*JJ * chx[i] );
        gs1x.push_back( p1x[i] - 1.*JJ * ch1x[i] );

        da.push_back( Dn[i+1] / Index + (double)(i+1) / SizeParam );
        db.push_back( Index * Dn[i+1] + (double)(i+1) / SizeParam );
        an[i] = (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) ;
        bn[i] = (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) ;
    }

}


void
SPHERE::LowFreqAnBn(complex128* an, complex128* bn)
{

  double LL, m2, x3, x4, x5, x6;

  m2          = Index * Index;
  LL          = (m2 - 1) / (m2 + 2);
  x3          = SizeParam * SizeParam * SizeParam;
  x4          = x3 * SizeParam;
  x5          = x4 * SizeParam;
  x6          = x5 * SizeParam;

  an[0] = (-2.*JJ * x3 / 3.) * LL - (2.*JJ * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  an[1] = (-1.*JJ * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  bn[0] = (-1.*JJ * x5 / 45.) * (m2 - 1.);
  bn[1] = 0. + 0.*JJ;

}

std::tuple<Cndarray,Cndarray>
SPHERE::S1S2(const ndarray Phi)
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
SPHERE::sFields(ndarray& Phi, ndarray& Theta, double R)
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

  complex128   propagator = E0 / (k * R) * exp(-JJ*k*R);

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, Polarization);

  std::tie(S1, S2) = this->S1S2(Phi);

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

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(EPhi, ETheta)  ;

}

std::tuple<Cndarray,Cndarray>
SPHERE::uFields(ndarray& Phi, ndarray& Theta, double R)
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

  complex128   propagator = E0 / (k * R) * exp(-JJ*k*R);

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, Polarization);

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
SPHERE::sS1S2(ndarray& Phi, ndarray& Theta)
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

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, Polarization);

  std::tie(S1, S2) = this->S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

  Structured(ThetaLength, PhiLength, S2Ptr, CosTerm, E0, EPhiPtr);

  Structured(ThetaLength, PhiLength, S1Ptr, SinTerm, E0, EThetaPtr);

  EPhi.resize({PhiLength,ThetaLength});
  ETheta.resize({PhiLength,ThetaLength});

  EPhi   = EPhi.attr("transpose")();
  ETheta = ETheta.attr("transpose")();

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(EPhi, ETheta)  ;
}


std::tuple<Cndarray,Cndarray>
SPHERE::uS1S2(ndarray& Phi, ndarray& Theta)
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

  PolarizationTerm(ThetaLength, ThetaPtr, CosTerm, SinTerm, Polarization);

  std::tie(S1, S2) = this->S1S2(Phi);

  complex128 * EPhiPtr    = (complex128*) ETheta.request().ptr,
             * EThetaPtr  = (complex128*) EPhi.request().ptr,
             * S1Ptr      = (complex128*) S1.request().ptr,
             * S2Ptr      = (complex128*) S2.request().ptr;

  Unstructured(ThetaLength, PhiLength, S2Ptr, CosTerm, E0, EPhiPtr);

  Unstructured(ThetaLength, PhiLength, S1Ptr, SinTerm, E0, EThetaPtr);

  free(CosTerm);
  free(SinTerm);
  return std::make_tuple(EPhi, ETheta)  ;

}



Cndarray
SPHERE::Dn(uint MaxOrder)
{
  double mx          = Index * SizeParam,
         M           = Index/nMedium;

  Cndarray Dn        = Cndarray(MaxOrder);

  complex128 *DnPtr  = (complex128*) Dn.request().ptr,
              numerator,
              denominator;

  for (uint order = 1; order < MaxOrder+1; order++)
  {
    numerator        = Mu * M*M * ( Xi(order, SizeParam) * Psi_p(order, SizeParam) - Xi_p(order, SizeParam) * Psi(order, SizeParam) );
    denominator      = Mu * M * Xi(order, SizeParam) * Psi_p(order, mx) - MuScat * Xi_p(order, SizeParam) * Psi(order, mx);
    DnPtr[order-1]   = numerator / denominator;
  }
  return Dn;
}


Cndarray
SPHERE::Cn(uint MaxOrder)
{
  double mx          = Index * SizeParam,
         M           = Index/nMedium;

  Cndarray Cn        = Cndarray(MaxOrder);

  complex128 *CnPtr  = (complex128*) Cn.request().ptr,
              numerator,
              denominator;

  for (uint order = 1; order < MaxOrder+1; order++)
  {
    numerator        = M * MuScat * ( Xi(order, SizeParam) * Psi_p(order, SizeParam) - Xi_p(order, SizeParam) * Psi(order, SizeParam) );
    denominator      = MuScat * Xi(order, SizeParam) * Psi_p(order, mx) - Mu * M * Xi_p(order, SizeParam) * Psi(order, mx);
    CnPtr[order-1]   = numerator / denominator;
  }
  return Cn;
}


Cndarray
SPHERE::Bn(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return bn;
}


Cndarray
SPHERE::An(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return an;
}


std::tuple<double, double, double, double, double, double, double>
SPHERE::GetEfficiencies()
{
    uint MaxOrder   = GetMaxOrder(SizeParam);

    complex128 * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
               * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128));

    this->ComputeAnBn(an, bn, MaxOrder);

    double Qsca   = GetQsca(an, bn, MaxOrder, SizeParam);
    double Qext   = GetQext(an, bn, MaxOrder, SizeParam);
    double g      = Getg(an, bn, MaxOrder, SizeParam, Qsca);
    double Qabs   = Qext - Qsca;
    double Qback  = GetQback(an, bn, MaxOrder, SizeParam);
    double Qpr    = Qext - g * Qsca;
    double Qratio = Qback / Qsca;

    free(an);
    free(bn);
    return std::make_tuple(Qsca, Qext, Qabs, Qback, Qratio, g, Qpr);
}






















// -
