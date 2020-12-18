#include "MieS1S2.hpp"
#include <vector>
#include <complex>
#include <boost/math/special_functions.hpp>
#include <cmath>
#include "Math.cpp"

#if __has_include("python3.8/Python.h")
#include "python3.8/Python.h"
#elif _has_include("python3.6/Python.h")
#include "python3.6/Python.h"
#endif


#define PI 3.14159265
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<std::vector<complex128>> iMatrix;



//REF: PhD Thesis   ON LIGHT SCATTERING BY NANOPARTICLES WITH CONVENTIONAL AND NON-CONVENTIONAL OPTICAL PROPERTIES
//REF: PhD Thesis https://www.google.com/url?sa=i&source=web&cd=&ved=2ahUKEwjvg4yF3cbtAhUPac0KHQj_BZkQ3YkBegQIARAE&url=http%3A%2F%2Frepositorio.unican.es%2Fxmlui%2Fbitstream%2Fhandle%2F10902%2F1566%2F2de8.BGCparteIcap2.pdf%3Fsequence%3D3&psig=AOvVaw18kz43dplVLIwhnDBQTTYI&ust=1607798286433047

static void
LowFrequencyMie_ab(const double m,
                   const double x,
                   iVec *an,
                   iVec *bn)
{
  const std::complex<double> j (0., 1.0);

  double LL, m2, x3, x4, x5, x6;
  complex128 a1, a2, b1, b2;

  m2 = m * m;
  LL = (m2 - 1) / (m2 + 2);
  x3 = x * x * x;
  x4 = x3 * x;
  x5 = x4 * x;
  x6 = x5 * x;

  a1 = (-2.*j * x3 / 3.) * LL - (2.*j * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
  a2 = (-1.*j * x5 / 15.) * (m2 - 1.) / (2. * m2 + 3.);
  b1 = (-1.*j * x5 / 45.) * (m2 - 1.);
  b2 = 0. + 0.*j;

  an->push_back(a1);
  an->push_back(a2);
  bn->push_back(b1);
  bn->push_back(b2);
}



static void
HighFrequencyMie_ab(const double m,
                    const double x,
                    const long unsigned int nmax,
                    const std::vector<double>* n,
                    iVec *an,
                    iVec *bn)

{
  const double mx = m * x;
  const double temp  = sqrt(0.5 * PI * x);
  const long unsigned int nmx = (long unsigned int) ( std::max( nmax, (long unsigned int) abs(mx) ) + 16 );
  iVec gsx, gs1x;
  iVec px, chx, p1x, ch1x, D, da, db;
  std::vector<double> Dn = std::vector<double>(nmx);
  std::complex<double> j (0., 1.0);

  p1x.push_back( sin(x) );
  ch1x.push_back( cos(x) );

  for (double i = nmx - 1; i > 1; i--)
  {
      Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );
  }

  for (long unsigned int i = 0; i < nmax; i++)
  {
    px.push_back(  temp * boost::math::cyl_bessel_j( (*n)[i] + 0.5, x ) );         //jv
    chx.push_back(-temp * boost::math::cyl_neumann(  (*n)[i] + 0.5, x ) );          //yv

    p1x.push_back(px[i]);
    ch1x.push_back(chx[i]);

    gsx.push_back( px[i] - 1.*j * chx[i] );
    gs1x.push_back( p1x[i] - 1.*j * ch1x[i] );

    D.push_back(Dn[i+1]);

    da.push_back( D[i] / m + (*n)[i] / x );
    db.push_back( m * D[i] + (*n)[i] / x );

    an->push_back( (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) );
    bn->push_back( (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) );
  }
}




static void
MiePiTau(const double mu,
         const long unsigned int nmax,
         iVec *pin,
         iVec *taun )

{
  (*pin)[0] = 1.;
  (*pin)[1] = 3. * mu;

  (*taun)[0] = mu;
  (*taun)[1] = 3.0 * cos(2. * acos(mu) );
  double temp = 0;
  for (long unsigned int i = 2; i < nmax; i++)
      {
       temp = (double)i;
       (*pin)[i] = ( (2. * temp + 1.) * ( mu * (*pin)[i-1] ) - (temp + 1.) * (*pin)[i-2] ) / temp;
       (*taun)[i] = (temp + 1.) * mu * (*pin)[i] - (temp + 2.) * (*pin)[i-1];
     }
}





static void
C_GetS1S2(const double m,
          const double x,
          const double*  phi,
          const long unsigned int lenght,
          complex128* S1S2)

{
    iVec *an = new iVec;
    iVec *bn = new iVec;

    std::vector<double> *n, *n2;

    const long unsigned int nmax = (int) round(2. + x + 4. * pow(x, 1./3.) );

    std::tie(n, n2) = Arrange(1, nmax + 1);

    (x < 0.5) ? LowFrequencyMie_ab(m, x, an, bn) : HighFrequencyMie_ab(m, x, nmax, n, an, bn);

    const long unsigned int anLength = an->size();

    iVec S1 = iVec(lenght) ;
    iVec S2 = iVec(lenght) ;

    iVec *pin = new iVec(nmax);
    iVec *taun = new iVec(nmax);
    complex128 j (0., 1.0);

    complex128 *temp0 = &S1S2[0], *temp1 = &S1S2[lenght] ;

    for (long unsigned int i = 0; i < lenght; i++){

        MiePiTau(cos( phi[i] ), nmax, pin, taun);

        for (long unsigned int k = 0; k < anLength ; k++){

            *temp0 += (*n2)[k] * ( (*an)[k] * (*pin)[k] +  (*bn)[k] * (*taun)[k] );
            *temp1 += (*n2)[k] * ( (*an)[k] * (*taun)[k] + (*bn)[k] * (*pin)[k] ) ;

          }
    temp0++ ;
    temp1++ ;
    }

    return;

}



static void
C_GetFields(const double m,
            const double x,
            const double* ThetaVec,
            const int Thetalenght,
            const double* PhiVec,
            const int Philenght,
            complex128* Parallel,
            complex128* Perpendicular,
            double Polarization
          )
{
  complex128* S1S2 = (complex128*) calloc(2 * Philenght , sizeof(complex128));

  const std::complex<double> j (0., 1.0) ;

  double temp0 ;
  complex128 temp2;

  C_GetS1S2(m, x, PhiVec, Philenght, S1S2) ;

  for (long unsigned int k=0; k < Thetalenght; k++ )
  {
    temp0 = *ThetaVec++ ;
    for (long unsigned int i=0; i < Thetalenght; i++ ){
      *Parallel++          = S1S2[i] * abs(cos(temp0 + Polarization)) ;
      *Perpendicular++     = S1S2[i + Philenght] * abs(sin(temp0 + Polarization)) ;
    }
  }

  free(S1S2) ;
  return;
}




static void
_C_GetFields(const double m,
            const double x,
            double* Theta,
            const double*  phi,
            const double* phiOnly,
            const long unsigned int Thetalenght,
            const long unsigned int Philenght,
            complex128* Parallel,
            complex128* Perpendicular,
            double Polarization
          )
{
  complex128* S1S2 = (complex128*) calloc(2 * Philenght , sizeof(complex128));

  const std::complex<double> j (0., 1.0) ;

  double temp0 ;
  complex128 temp2;

  C_GetS1S2(m, x, phiOnly, Philenght, S1S2) ;
/*
  long unsigned int i = 0;
  for (long unsigned int k=0; k < Thetalenght; k++ )
  {
    temp0 = *Theta++ ;

    if (k % Philenght == 0) i++;

    *Parallel++          = S1S2[i] * abs(cos(temp0 + Polarization)) ;
    *Perpendicular++     = S1S2[i + Philenght] * abs(sin(temp0 + Polarization)) ;



  }

  free(S1S2) ;
  */
  return;
}





static void
C_GetFieldsNoPolarization(const double m,
                          const double x,
                          const double* ThetaVec,
                          const int Thetalenght,
                          const double* PhiVec,
                          const int Philenght,
                          complex128* Parallel,
                          complex128* Perpendicular)
{
  complex128* S1S2 = (complex128*) calloc(2 * Philenght , sizeof(complex128));

  const std::complex<double> j (0., 1.0) ;

  double temp0 ;
  double temp1 = 1./sqrt(2.) ;
  complex128 temp2;

  C_GetS1S2(m, x, PhiVec, Philenght, S1S2) ;

  for (long unsigned int k=0; k < Thetalenght; k++ )
  {
    temp0 = *ThetaVec++ ;
    for (long unsigned int i=0; i < Thetalenght; i++ ){
      *Parallel++          = S1S2[i] * temp1 ;
      *Perpendicular++     = S1S2[i + Philenght] * temp1 ;
    }
  }

  free(S1S2) ;
  return;
}





// -
