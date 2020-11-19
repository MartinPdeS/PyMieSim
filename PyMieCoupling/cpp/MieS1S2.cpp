//#include "MieS1S2.h"
#include <vector>
#include <complex>
#include<iostream>
#include <tuple>
#include <algorithm>
#include <boost/math/special_functions.hpp> // For gamma function.
#include <boost/math/special_functions/bessel.hpp>
#include <math.h>
#include <cmath>
#include "Math.cpp"
#define PI 3.14159265

typedef std::vector<std::complex<double>> iVec;
typedef std::complex<double> complex128;







std::tuple<std::vector<complex128> , std::vector<complex128>> LowFrequencyMie_ab(double m,
                                                                                 double x
                                                                                )
{
  std::complex<double> j (0., 1.0);

  double LL, m2, x3, x4, x5, x6;
  complex128 a1, a2, b1, b2;
  std::vector<complex128> an, bn;

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

  an.push_back(a1); an.push_back(a2);
  bn.push_back(b1); bn.push_back(b2);


  return std::make_tuple(an, bn);
}



std::tuple<std::vector<complex128> , std::vector<complex128>> HighFrequencyMie_ab(double m,
                                                                                  double x,
                                                                                  int nmax,
                                                                                  std::vector<double> n
                                                                                  )

{
  double mx = m * x;
  double temp  = sqrt(0.5 * PI * x);
  int nmx = (int) ( std::max( nmax, (int) abs(mx) ) + 16 );
  std::vector<complex128> an, bn, gsx, gs1x;
  std::vector<complex128> px, chx, p1x, ch1x, D, da, db;
  std::vector<double> Dn = std::vector<double>(nmx);
  std::complex<double> j (0., 1.0);


  p1x.push_back( sin(x) );
  ch1x.push_back( cos(x) );

  for (double i = nmx - 1; i < 1; i--)
  {
      Dn[i-1] = (i / mx) - ( 1 / (Dn[i] + i/mx) );
  }

  for (int i = 0; i < nmax; i++)
  {
    px.push_back(  temp * boost::math::cyl_bessel_j( n[i] + 0.5, x ) );         //jv
    chx.push_back(-temp * boost::math::cyl_neumann( n[i] + 0.5, x ) );          //yv

    gsx.push_back( px[i] - 1.*j * chx[i] );
    gs1x.push_back( p1x[i] - 1.*j * ch1x[i] );

    D.push_back(Dn[i+1]);

    da.push_back( D[i] / m + n[i] / x );
    db.push_back( m * D[i] + n[i] / x );

    an.push_back( (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) );
    bn.push_back( (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) );

  }


  return std::make_tuple(an, bn);

}






std::tuple<std::vector<complex128> , std::vector<complex128>> MiePiTau(double m,
                                                                       double mu,
                                                                       double x,
                                                                       int nmax)

{
  std::vector<complex128> pin, taun;
  pin.push_back(1.); pin.push_back(3. * mu);
  taun.push_back(mu); taun.push_back(3.0 * cos(2. * acos(mu) ) );

  for (int i = 0; i < nmax; i++)
  {
    pin.push_back( ( (2 * (double)i + 1) * ( mu * pin[i-1] ) - ((double)i + 1) * pin[i-2] ) / (double)i );
    taun.push_back( ((double)i + 1) * mu * pin[i] - ((double)i + 2) * pin[i-1] );
  }

  return std::make_tuple(pin, taun);

}









std::tuple<std::vector<complex128> , std::vector<complex128>> test(double m, double x, std::vector<double> phi)
{
    std::vector<complex128> an, bn, SS1, SS2, pin, taun, S1, S2;
    std::vector<double> n, n2;

    int nmax = (int) (2. + x + 4. * pow(x, 1./3.) );
    double mu = 0.;

    std::tie(an, bn) = LowFrequencyMie_ab(m, x);

    std::tie(n, n2) = Arrange(1, nmax + 1);

    for (int i = 0; i < phi.size(); i++){
        mu = phi[i];
        std::tie(pin, taun) = MiePiTau(m, mu, x, nmax);

        for (int i = 0; i < an.size(); i++){
            SS1.push_back( n2[i] * ( an[i] * pin[i] + bn[i] * taun[i] ) );
            SS2.push_back( n2[i] * ( an[i] * taun[i] + bn[i] * pin[i] ) );
            }

        S1.push_back(Sum(SS1)); S2.push_back(Sum(SS2));

    }

  return std::make_tuple(S1, S2);

}












// -
