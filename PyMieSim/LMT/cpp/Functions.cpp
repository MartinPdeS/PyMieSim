#include <vector>
#include <complex>
#include <cmath>


typedef std::complex<double> complex128;

#define J std::complex<double>(0.0,1.0)

#define PI (double)3.14159265358979323846264338





void
LoopUnstructured(int         PhiLength,
                 int         ThetaLength,
                 complex128  propagator,
                 double     *PhiPtr,
                 double     *ThetaPtr,
                 complex128 *EPhiPtr,
                 complex128 *EThetaPtr,
                 complex128 *s1s2,
                 bool        Polarized,
                 double      Polarization)

{
  int w;
  double temp0;

  for (auto p=0; p < PhiLength; p++ )
  {

        temp0 = ThetaPtr[p] ;
        EPhiPtr[w] = J* propagator * s1s2[p] * (complex128) abs(cos(temp0 + Polarization));
        EThetaPtr[w] =- propagator * s1s2[p + PhiLength] * (complex128) abs(sin(temp0 + Polarization));
        w++;

  }

}



void
LoopStructured(int         PhiLength,
               int         ThetaLength,
               complex128  propagator,
               double     *PhiPtr,
               double     *ThetaPtr,
               complex128 *EPhiPtr,
               complex128 *EThetaPtr,
               complex128 *s1s2,
               bool        Polarized,
               double      Polarization)

{
  int w = 0;
  double temp0;

  for (auto p=0; p < PhiLength; p++ )
  {
    for (auto t=0; t < ThetaLength; t++ )
        {
          temp0 = ThetaPtr[t] ;
          EPhiPtr[w] = J* propagator * s1s2[p] * (complex128) abs(cos(temp0 + Polarization));
          EThetaPtr[w] =- propagator * s1s2[p + PhiLength] * (complex128) abs(sin(temp0 + Polarization));
          w++;
       }
  }

}
