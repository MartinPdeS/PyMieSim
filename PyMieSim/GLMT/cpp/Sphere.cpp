#include <iostream>


class _SPHERE{

private:
  bool          Polarized;

  double        Diameter,
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
                GetQext();

    void        ComputeAnBn(complex128* an, complex128* bn, uint MaxOrder),
                LowFreqAnBn(complex128* an, complex128* bn),
                HighFreqAnBn(complex128* an, complex128* bn, uint MaxOrder),
                IsPolarized(),
                ComputePrefactor(double* prefactor, uint MaxOrder),
                PolarizationTerm(uint ThetaLength, double * ThetaPtr, double * CosTerm, double * SinTerm);

    inline void LoopBSC(Cndarray& BSC, double* prefactor, complex128* an, complex128* bn, complex128* pin, complex128* taun, complex128& s1, complex128& s2, double& Theta),
                MiePiTau(double mu, uint MaxOrder, complex128 *pin, complex128 *taun);

    Cndarray    BSC;


    public:

        std::tuple<double, double, double> GetEfficiencies();
        Cndarray                           An(uint MaxOrder);
        Cndarray                           Bn(uint MaxOrder);
        Cndarray                           Cn(uint MaxOrder);
        Cndarray                           Dn(uint MaxOrder);
        std::tuple<Cndarray,Cndarray>      sS1S2(ndarray& Phi, ndarray& Theta);
        std::tuple<Cndarray,Cndarray>      uS1S2(ndarray& Phi, ndarray& Theta);
        std::tuple<Cndarray,Cndarray>      sFields(ndarray& Phi, ndarray& Theta, double R);
        std::tuple<Cndarray,Cndarray>      uFields(ndarray& Phi, ndarray& Theta, double R);


  _SPHERE(double Index,
         double Diameter,
         double Wavelength,
         double nMedium,
         double Polarization,
         double E0,
         Cndarray BSC)
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
          this->BSC           = BSC;
          this->IsPolarized();
        }

        ~_SPHERE(){  }

};


void
_SPHERE::PolarizationTerm(uint ThetaLength, double * ThetaPtr, double * CosTerm, double * SinTerm)
{
  if (this->Polarized==true)
  {
    for (uint t = 0; t < ThetaLength; t++)
    {
        CosTerm[t] = abs(cos(Polarization + ThetaPtr[t])) ;
        SinTerm[t] = abs(sin(Polarization + ThetaPtr[t])) ;
    }
  }
  else
  {
    const double term = 1./sqrt(2);
    for (uint t = 0; t < ThetaLength; t++)
    {
        CosTerm[t] = term ;
        SinTerm[t] = term ;
    }
  }
}


void
_SPHERE::IsPolarized()
{
  if (Polarization==-1.){Polarized=false;}
  else                  {Polarized=true;}
}



std::tuple<Cndarray,Cndarray>
_SPHERE::sS1S2(ndarray& Phi, ndarray& Theta)
{

  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0],
               BSCLength    = BSC.request().shape[0],
               MaxOrder     = (uint)((complex128 *)BSC.request().ptr)[BSCLength-1].real();

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * PhiPtr       = (double*) Phi.request().ptr,
             * prefactor    = (double*) calloc(MaxOrder+1, sizeof(double));

  Cndarray     S1       = Cndarray(PhiLength*ThetaLength),
               S2         = Cndarray(PhiLength*ThetaLength);

  complex128 * an           = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * bn           = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * pin          = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * taun         = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * S1Ptr        = (complex128*) S1.request().ptr,
             * S2Ptr        = (complex128*) S2.request().ptr,
               s1           = 0.,
               s2           = 0.;

  this->ComputeAnBn(an, bn, MaxOrder);

  this->ComputePrefactor(prefactor, MaxOrder+1);

   uint index     = 0;

   for(uint p = 0; p < PhiLength; p++)
     {
       for(uint t = 0; t < ThetaLength; t++)
         {
            this->MiePiTau( cos( PhiPtr[p]-PI/2 ), MaxOrder+1, pin, taun );

            this->LoopBSC(BSC, prefactor, an, bn, pin, taun, s1, s2, ThetaPtr[t]);

            S1Ptr[index]   = s1;
            S2Ptr[index]   = s2;
            index++;
         }
    }

  S1.resize({PhiLength,ThetaLength});
  S2.resize({PhiLength,ThetaLength});

  free(an);
  free(bn);
  free(pin);
  free(taun);
  free(prefactor);

  return std::make_tuple(S1, S2);
}


std::tuple<Cndarray,Cndarray>
_SPHERE::uS1S2(ndarray& Phi, ndarray& Theta)
{

  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0],
               BSCLength    = BSC.request().shape[0],
               MaxOrder     = (uint)((complex128 *)BSC.request().ptr)[BSCLength-1].real();

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * PhiPtr       = (double*) Phi.request().ptr,
             * prefactor    = (double*) calloc(MaxOrder+1, sizeof(double));

  Cndarray     S1       = Cndarray(PhiLength*ThetaLength),
               S2         = Cndarray(PhiLength*ThetaLength);

  complex128 * an           = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * bn           = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * pin          = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * taun         = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * S1Ptr        = (complex128*) S1.request().ptr,
             * S2Ptr        = (complex128*) S2.request().ptr,
               s1           = 0.,
               s2           = 0.;

  this->ComputeAnBn(an, bn, MaxOrder);

  this->ComputePrefactor(prefactor, MaxOrder+1);

   for(uint p = 0; p < PhiLength; p++)
     {
        this->MiePiTau( cos( PhiPtr[p]-PI/2 ), MaxOrder+1, pin, taun );
        this->LoopBSC(BSC, prefactor, an, bn, pin, taun, s1, s2, ThetaPtr[p]);

        S1Ptr[p]   = s1;
        S2Ptr[p]   = s2;
     }

  S1.resize({PhiLength,ThetaLength});
  S2.resize({PhiLength,ThetaLength});

  free(an);
  free(bn);
  free(pin);
  free(taun);
  free(prefactor);
  return std::make_tuple(S1, S2);
}


std::tuple<Cndarray,Cndarray>
_SPHERE::sFields(ndarray& Phi, ndarray& Theta, double R)
{
  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0],
               index        = 0;

  Cndarray     ETheta, EPhi;

  std::tie(EPhi, ETheta) = this->sS1S2(Phi, Theta);

  complex128 * EPhiPtr      = (complex128*) ETheta.request().ptr,
             * EThetaPtr    = (complex128*) EPhi.request().ptr,
               propagator   = E0 / (k * R) * exp(-JJ*k*R);

  for(uint p = 0; p < PhiLength; p++)
    {
      for(uint t = 0; t < ThetaLength; t++)
        {
           EPhiPtr[index]   *= JJ  * propagator;
           EThetaPtr[index] *= - propagator;
           index++;
        }
   }
  return std::make_tuple(ETheta, EPhi);
}


std::tuple<Cndarray,Cndarray>
_SPHERE::uFields(ndarray& Phi, ndarray& Theta, double R)
{
  uint         PhiLength    = Phi.request().shape[0];

  Cndarray     ETheta, EPhi;

  std::tie(EPhi, ETheta) = this->uS1S2(Phi, Theta);

  complex128 * EPhiPtr      = (complex128*) ETheta.request().ptr,
             * EThetaPtr    = (complex128*) EPhi.request().ptr,
               propagator   = E0 / (k * R) * exp(-JJ*k*R);

  for(uint p = 0; p < PhiLength; p++)
    {
       EPhiPtr[p]   *= JJ * propagator;
       EThetaPtr[p] *= -  propagator;
    }
  return std::make_tuple(ETheta, EPhi);
}


inline void
_SPHERE::LoopBSC(Cndarray& BSC,
                double* prefactor,
                complex128* an,
                complex128* bn,
                complex128* pin,
                complex128* taun,
                complex128& s1,
                complex128& s2,
                double&     Theta)
{
  uint BSCLength      = BSC.request().shape[0];

  complex128 * BSCPtr =  (complex128*)BSC.request().ptr;

  for (uint b = 0; b < BSCLength; b++)
    {
        int n         = (int)BSCPtr[b].real();
        int m         = (int)BSCPtr[b + BSCLength].real();
        complex128 TE = BSCPtr[b + BSCLength*2];
        complex128 TM = BSCPtr[b + BSCLength*3];

        complex128 _exp   = exp( JJ * (double)m * Theta + Polarization );

        s1 += prefactor[n]*(JJ * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;
        s2 += prefactor[n]*(JJ * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;
   }
}


double
_SPHERE::GetQsca()
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
_SPHERE::GetQext()
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
_SPHERE::GetEfficiencies()
{
    double Qsca = GetQsca();
    double Qext = GetQext();
    double Qabs = Qext - Qsca;

    return std::make_tuple(Qsca, Qext, Qabs);
}


void
_SPHERE::ComputeAnBn(complex128* an, complex128* bn, uint MaxOrder)
{
  if (SizeParam < 0.5){LowFreqAnBn(an, bn) ; }
  else                {HighFreqAnBn(an, bn, MaxOrder) ; }
}


void
_SPHERE::HighFreqAnBn(complex128* an, complex128* bn, uint MaxOrder)
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
_SPHERE::LowFreqAnBn(complex128* an, complex128* bn)
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


Cndarray
_SPHERE::Dn(uint MaxOrder)
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
_SPHERE::Cn(uint MaxOrder)
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
_SPHERE::Bn(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return bn;
}


Cndarray
_SPHERE::An(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return an;
}


inline void
_SPHERE::MiePiTau(double        mu,
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
_SPHERE::ComputePrefactor(double* prefactor, uint MaxOrder)
{
   for (uint m = 0; m < MaxOrder ; m++)
   {
      prefactor[m] = (double) ( 2 * (m+1) + 1 ) / ( (m+1) * ( (m+1) + 1 ) );
   }
}











































//-
