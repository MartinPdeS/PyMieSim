#include <iostream>


class BASE{
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

public:

    void  ComputePrefactor(double* prefactor, uint MaxOrder),
          IsPolarized(),
          Structured(uint ThetaLength, uint PhiLength, complex128 *array0, double *array1, complex128  scalar, complex128 *output),
          Unstructured(uint ThetaLength, uint PhiLength, complex128 *array0, double *array1, complex128  scalar, complex128 *output),
          ComputeAnBn(complex128* an, complex128* bn, uint MaxOrder),
          PolarizationTerm(uint ThetaLength, double * ThetaPtr, double * CosTerm, double * SinTerm);

    inline void MiePiTau(double mu, uint MaxOrder, complex128 *pin, complex128 *taun);


    BASE(){}

    virtual ~BASE(){ }

};


void
BASE::PolarizationTerm(uint ThetaLength, double * ThetaPtr, double * CosTerm, double * SinTerm)
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
BASE::IsPolarized()
{
  if (Polarization==-1.){Polarized=true;}
  else                  {Polarized=true;}
}



void
BASE::Unstructured(uint        ThetaLength,
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


void
BASE::Structured(uint        ThetaLength,
                 uint        PhiLength,
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


inline void
BASE::MiePiTau(double      mu,
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
BASE::ComputePrefactor(double* prefactor, uint MaxOrder)
{
   for (uint m = 0; m < MaxOrder ; m++)
   {
      prefactor[m] = (double) ( 2 * (m+1) + 1 ) / ( (m+1) * ( (m+1) + 1 ) );
   }
}








// -
