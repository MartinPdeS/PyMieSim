#include <iostream>


class _CYLINDER{

private:
  bool          Polarized;

  uint          BSCLength,
                MaxOrder;

  complex128  * BSCPtr;

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

    void        ComputeAnBn(complex128* an, complex128* bn),
                LowFreqAnBn(complex128* an, complex128* bn),
                HighFreqAnBn(complex128* an, complex128* bn, uint MaxOrder),
                ComputeExpTerm(complex128 * expTerm, int ThetaLength, double* ThetaPtr),
                ComputeMuTerm(double * muTerm, int PhiLength, double* PhiPtr),
                ComputeTerms(complex128 * Term1, complex128 * Term2, complex128 * an, complex128 * bn, double * prefactor),
                PolarizationTerm(uint ThetaLength, double * ThetaPtr, double * CosTerm, double * SinTerm),
                ComputePrefactor(double* prefactor),
                IsPolarized();

    inline void LoopBSC(double* prefactor, complex128* an, complex128* bn, double* pin, double* taun, complex128& s1, complex128& s2, double& Theta),
                MiePiTau(double &mu, double *pin, double *taun);

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


  _CYLINDER(double   Index,
          double   Diameter,
          double   Wavelength,
          double   nMedium,
          double   Polarization,
          double   E0,
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
          this->BSCLength     = BSC.request().shape[0];
          this->BSCPtr        = (complex128*) BSC.request().ptr;
          this->MaxOrder      = (uint)((complex128 *)BSC.request().ptr)[BSCLength-1].real();
          this->IsPolarized();
        }

        ~_CYLINDER(){  }

};


void
_CYLINDER::PolarizationTerm(uint ThetaLength, double * ThetaPtr, double * CosTerm, double * SinTerm)
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
_CYLINDER::IsPolarized()
{
  if (Polarization==-1.){Polarized=false;}
  else                  {Polarized=true;}
}



std::tuple<Cndarray,Cndarray>
_CYLINDER::sS1S2(ndarray& Phi, ndarray& Theta)
{

  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0],
               n            = 0;

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * PhiPtr       = (double*) Phi.request().ptr,
             * prefactor    = (double*) calloc(MaxOrder, sizeof(double)),
             * muTerm       = (double*) calloc(ThetaLength, sizeof(double)),
             * pin          = (double*) calloc(MaxOrder, sizeof(double)),
             * taun         = (double*) calloc(MaxOrder, sizeof(double)),
               m            = 0.;

  Cndarray     S1           = Cndarray(PhiLength*ThetaLength),
               S2           = Cndarray(PhiLength*ThetaLength);

  complex128 * an           = (complex128*) calloc(MaxOrder, sizeof(complex128)),
             * bn           = (complex128*) calloc(MaxOrder, sizeof(complex128)),
             * S1Ptr        = (complex128*) S1.request().ptr,
             * S2Ptr        = (complex128*) S2.request().ptr,
               s1           = 0.,
               s2           = 0.,
             * expTerm      = (complex128*) calloc(BSCLength*ThetaLength, sizeof(complex128)),
             * Term1        = (complex128*) calloc(BSCLength, sizeof(complex128)),
             * Term2        = (complex128*) calloc(BSCLength, sizeof(complex128)),
               term3;

   this->ComputeAnBn(an, bn);

   this->ComputePrefactor(prefactor);

   this->ComputeMuTerm(muTerm, PhiLength, PhiPtr);

   this->ComputeExpTerm(expTerm, ThetaLength, ThetaPtr);

   this->ComputeTerms(Term1, Term2, an, bn, prefactor);

   for(uint p = 0; p < PhiLength; p++)
     {
       this->MiePiTau( muTerm[p], pin, taun );

       for(uint t = 0; t < ThetaLength; t++)
         {
           s1=0.;s2=0.;
           for (uint b = 0; b < BSCLength; b++)
             {
                 n             = (uint)BSCPtr[b].real();
                 m             = BSCPtr[b + BSCLength].real();
                 term3         = expTerm[t*BSCLength + b];

                 s1           += (Term1[b] * taun[n] + m * Term2[b] * pin[n] )  * term3;
                 s2           += (Term1[b] * m * pin[n] +  Term2[b] * taun[n] ) * term3;
            }
            *S1Ptr   = s1;
            *S2Ptr   = s2;

             S1Ptr++; S2Ptr++;


         }
    }

  S1.resize({PhiLength,ThetaLength});
  S2.resize({PhiLength,ThetaLength});

  free(an);
  free(bn);
  free(pin);
  free(taun);
  free(prefactor);
  free(expTerm);
  free(muTerm);
  free(Term1);
  free(Term2);

  return std::make_tuple(S1, S2);
}

void
_CYLINDER::ComputeTerms(complex128 * Term1, complex128 * Term2, complex128 * an, complex128 * bn, double * prefactor)
{
    uint n;
    for (uint b = 0; b < BSCLength; b++)
    {
        n             = (uint)BSCPtr[b].real();
        Term1[b]      = prefactor[n] * JJ * bn[n] * BSCPtr[b + BSCLength*2] ;
        Term2[b]      = prefactor[n] * an[n] * BSCPtr[b + BSCLength*3];
    }
}

void
_CYLINDER::ComputeMuTerm(double * muTerm, int PhiLength, double* PhiPtr)
{
    for(uint p = 0; p < PhiLength; p++)
    {
      muTerm[p] = cos( PhiPtr[p]-PI/2 ) ;
    }
}


void
_CYLINDER::ComputeExpTerm(complex128 * expTerm, int ThetaLength, double* ThetaPtr)
{
  for(uint t = 0; t < ThetaLength; t++)
    {
      for (uint b = 0; b < BSCLength; b++)
        {
           double m       = BSCPtr[b + BSCLength].real();
           expTerm[t*BSCLength + b] = exp( JJ * m * ThetaPtr[t] + Polarization );
        }
    }
}


std::tuple<Cndarray,Cndarray>
_CYLINDER::uS1S2(ndarray& Phi, ndarray& Theta)
{

  uint         PhiLength    = Phi.request().shape[0],
               ThetaLength  = Theta.request().shape[0];

  double     * ThetaPtr     = (double*) Theta.request().ptr,
             * PhiPtr       = (double*) Phi.request().ptr,
             * prefactor    = (double*) calloc(MaxOrder, sizeof(double)),
             * muTerm       = (double*) calloc(ThetaLength, sizeof(double)),
             * pin          = (double*) calloc(MaxOrder, sizeof(double)),
             * taun         = (double*) calloc(MaxOrder, sizeof(double));

  Cndarray     S1       = Cndarray(PhiLength*ThetaLength),
               S2         = Cndarray(PhiLength*ThetaLength);

  complex128 * an           = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * bn           = (complex128*) calloc(MaxOrder+1, sizeof(complex128)),
             * S1Ptr        = (complex128*) S1.request().ptr,
             * S2Ptr        = (complex128*) S2.request().ptr,
               s1           = 0.,
               s2           = 0.;

  this->ComputeAnBn(an, bn);

  this->ComputePrefactor(prefactor);

  for(uint p = 0; p < PhiLength; p++){muTerm[p] = cos( PhiPtr[p]-PI/2 ) ;}

   for(uint p = 0; p < PhiLength; p++)
     {
        s1=0.;s2=0.;
        this->MiePiTau( muTerm[p], pin, taun );
        this->LoopBSC(prefactor, an, bn, pin, taun, s1, s2, ThetaPtr[p]);

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
_CYLINDER::sFields(ndarray& Phi, ndarray& Theta, double R)
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
           EThetaPtr[index]   *= JJ  * propagator;
           EPhiPtr[index]     *= -1. * propagator;
           index++;
        }
   }
  return std::make_tuple(EPhi, ETheta);
}


std::tuple<Cndarray,Cndarray>
_CYLINDER::uFields(ndarray& Phi, ndarray& Theta, double R)
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
  return std::make_tuple(EPhi, ETheta);
}


inline void
_CYLINDER::LoopBSC(double*     prefactor,
                 complex128 *an,
                 complex128 *bn,
                 double     *pin,
                 double     *taun,
                 complex128& s1,
                 complex128& s2,
                 double&     Theta)
{
  for (uint b = 0; b < BSCLength; b++)
    {
        int n         = (int)BSCPtr[b].real();
        double m      = (double)BSCPtr[b + BSCLength].real();
        complex128 TE = BSCPtr[b + BSCLength*2];
        complex128 TM = BSCPtr[b + BSCLength*3];

        complex128 _exp   = exp( JJ * m * Theta + Polarization );

        s1 += prefactor[n]*(JJ * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;
        s2 += prefactor[n]*(JJ * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;
   }
}


double
_CYLINDER::GetQsca()
{
    uint MaxOrder   = GetMaxOrder(SizeParam);

    complex128 * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
               * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128));

    this->ComputeAnBn(an, bn);

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
_CYLINDER::GetQext()
{
    uint MaxOrder   = GetMaxOrder(SizeParam);

    complex128 * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
               * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128));

    this->ComputeAnBn(an, bn);

    complex128 temp = 0.;

    for(uint it = 0; it < MaxOrder; ++it){ temp += ( 2.*(double)(it+1) + 1.) * ( std::real( an[it] + an[it] ) ); }

    free(an);
    free(bn);

    return 2. / (SizeParam * SizeParam) * std::real(temp);
}


std::tuple<double, double, double>
_CYLINDER::GetEfficiencies()
{
    double Qsca = GetQsca();
    double Qext = GetQext();
    double Qabs = Qext - Qsca;

    return std::make_tuple(Qsca, Qext, Qabs);
}


inline void
_CYLINDER::MiePiTau(double  &mu,
                    double  *pin,
                    double  *taun)

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
_CYLINDER::ComputePrefactor(double* prefactor)
{
   for (uint m = 0; m < MaxOrder ; m++)
   {
      prefactor[m] = (double) ( 2 * (m+1) + 1 ) / ( (m+1) * ( (m+1) + 1 ) );
   }
}



void
_CYLINDER::ComputeAnBn(complex128* anPtr, complex128* bnPtr)
{
  return HighFreqAnBn(anPtr, bnPtr, MaxOrder) ;
}


void
_CYLINDER::HighFreqAnBn(complex128* anPtr, complex128* bnPtr, uint MaxOrder)
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
_CYLINDER::LowFreqAnBn(complex128* anPtr, complex128* bnPtr){}


Cndarray
_CYLINDER::Bn(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return bn;
}


Cndarray
_CYLINDER::An(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return an;
}




//-
