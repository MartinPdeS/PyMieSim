//#include "MieS1S2.h"
#include <vector>
#include <complex>
#include <boost/math/special_functions.hpp>
#include <cmath>
#include "Math.cpp"
//#include "python3.6/Python.h"
#if __has_include("python3.8/Python.h")
#include "python3.8/Python.h"
#elif _has_include("python3.6/Python.h")
#include "python3.8/Python.h"
#endif


//#include "/usr/lib/python3/dist-packages/numpy/core/include/numpy/arrayobject.h"
//#include "numpy/arrayobject.h"
//#include "list_of_ndarrays_lib.h"

#define PI 3.14159265
typedef std::vector<std::complex<double>> iVec;
typedef std::complex<double> complex128;





void LowFrequencyMie_ab(const double m,
                        const double x,
                        std::vector<complex128> *an,
                        std::vector<complex128> *bn)
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



void HighFrequencyMie_ab(const double m,
                         const double x,
                         const long unsigned int nmax,
                         const std::vector<double>* n,
                         std::vector<complex128> *an,
                         std::vector<complex128> *bn
                                            )

{
  const double mx = m * x;
  const double temp  = sqrt(0.5 * PI * x);
  const long unsigned int nmx = (long unsigned int) ( std::max( nmax, (long unsigned int) abs(mx) ) + 16 );
  std::vector<complex128> gsx, gs1x;
  std::vector<complex128> px, chx, p1x, ch1x, D, da, db;
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




void MiePiTau(const double mu,
              const long unsigned int nmax,
              std::vector<complex128> *pin,
              std::vector<complex128> *taun )

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





static std::pair<iVec, iVec> Cwrapper(const double m,
                                const double x,
                                const std::vector<double> phi)
{

    std::vector<complex128> *an = new std::vector<complex128>;
    std::vector<complex128> *bn = new std::vector<complex128>;

    std::vector<double> *n, *n2;

    const long unsigned int nmax = (int) round(2. + x + 4. * pow(x, 1./3.) );

    std::tie(n, n2) = Arrange(1, nmax + 1);

    if (x < 0.5){ LowFrequencyMie_ab(m, x, an, bn); }
    else{ HighFrequencyMie_ab(m, x, nmax, n, an, bn); }

    const long int anLength = an->size();
    const long int lenght = phi.size();

    std::vector<complex128> S1 = std::vector<complex128>(lenght) ;

    std::vector<complex128> S2 = std::vector<complex128>(lenght) ;

    std::vector<complex128> *pin = new std::vector<complex128>(nmax);
    std::vector<complex128> *taun = new std::vector<complex128>(nmax);
    std::complex<double> j (0., 1.0);

    PyObject* Py_S1 = PyList_New(0);


    for (long unsigned int i = 0; i < lenght; i++){

        MiePiTau(cos( phi[i] ), nmax, pin, taun);
        for (long unsigned int k = 0; k < anLength ; k++){


            S1[i] += (*n2)[k] * ( (*an)[k] * (*pin)[k] +  (*bn)[k] * (*taun)[k] );
            S2[i] += (*n2)[k] * ( (*an)[k] * (*taun)[k] + (*bn)[k] * (*pin)[k] ) ;

          }
    }

    return std::make_pair(S1, S2);


}



int main()
{

  return 1;
}









// -
