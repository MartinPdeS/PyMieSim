#include <vector>
#include <complex>
#include <boost/math/special_functions.hpp>
#include <cmath>
#include "Math.cpp"

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;


int GetMaxOrder(double SizeParam) {return (int) (2 + SizeParam + 4 * pow(SizeParam,1./3.)); }


std::tuple<iVec, iVec>
LowFrequencyMie_ab(const double m,
                   const double x)
{
  iVec an, bn ;
  const complex128 j (0., 1.0);

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

  an.push_back(a1);
  an.push_back(a2);
  bn.push_back(b1);
  bn.push_back(b2);
  return std::make_tuple(an, bn);
}



std::tuple<iVec, iVec>
HighFrequencyMie_ab(const double m,
                    const double x,
                    int          OrderMax,
                    const        Vec n)

{
  iVec an, bn ;
  const double mx = m * x,
               temp  = sqrt(0.5 * PI * x);
  const int nmx = (int) ( std::max( OrderMax, (int) abs(mx) ) + 16 );
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

  for (auto i = 0; i < OrderMax; i++)
  {
    px.push_back(  temp * boost::math::cyl_bessel_j( n[i] + 0.5, x ) );         //jv
    chx.push_back(-temp * boost::math::cyl_neumann(  n[i] + 0.5, x ) );          //yv

    p1x.push_back(px[i]);
    ch1x.push_back(chx[i]);

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




void
HighFrequencyMie_ab(const double               m,
                    const double               x,
                    const long unsigned int    OrderMax,
                    const Vec n,
                    iVec                      an,
                    iVec                      bn)

{
  const double mx = m * x;
  const double temp  = sqrt(0.5 * PI * x);
  const long unsigned int nmx = (long unsigned int) ( std::max( OrderMax, (long unsigned int) abs(mx) ) + 16 );
  iVec gsx, gs1x;
  iVec px, chx, p1x, ch1x, D, da, db;
  Vec Dn = Vec(nmx);
  complex128 j (0., 1.0);

  p1x.push_back( sin(x) );
  ch1x.push_back( cos(x) );

  for (double i = nmx - 1; i > 1; i--)
  {
      Dn[i-1] = (i / mx) - ( 1. / (Dn[i] + i/mx) );
  }

  for (long unsigned int i = 0; i < OrderMax; i++)
  {
    px.push_back(  temp * boost::math::cyl_bessel_j( n[i] + 0.5, x ) );         //jv
    chx.push_back(-temp * boost::math::cyl_neumann(  n[i] + 0.5, x ) );          //yv

    p1x.push_back(px[i]);
    ch1x.push_back(chx[i]);

    gsx.push_back( px[i] - 1.*j * chx[i] );
    gs1x.push_back( p1x[i] - 1.*j * ch1x[i] );

    D.push_back(Dn[i+1]);

    da.push_back( D[i] / m + n[i] / x );
    db.push_back( m * D[i] + n[i] / x );

    an.push_back( (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) );
    bn.push_back( (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) );
  }
}


std::tuple<iVec, iVec>
MiePiTau(const double            mu,
         const long unsigned int OrderMax)

{
  iVec pin = iVec(OrderMax), taun = iVec(OrderMax);
  pin[0] = 1.;
  pin[1] = 3. * mu;

  taun[0] = mu;
  taun[1] = 3.0 * cos(2. * acos(mu) );
  double n = 0;
  for (long unsigned int i = 2; i < OrderMax; i++)
      {
       n = (double)i;

       pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

       taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
     }
     return std::make_tuple(pin, taun);
}




void
MiePiTau(const double            mu,
         const long unsigned int OrderMax,
         iVec                    pin,
         iVec                    taun )

{
  pin[0] = 1.;
  pin[1] = 3. * mu;

  taun[0] = mu;
  taun[1] = 3.0 * cos(2. * acos(mu) );
  double n = 0;
  for (long unsigned int i = 2; i < OrderMax; i++)
      {
       n = (double)i;

       pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

       taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
     }
}


static double
C_Qsca(iVec           an,
       iVec           bn,
       const double   x,
       const Vec      n)

{
     double Qsca = 2. / (x * x);
     complex128 temp;
     for(long unsigned int it = 0; it < an.size(); ++it)
     {
       temp += (2.*n[it]+1) * (   std::real( an[it] ) * std::real( an[it] )
                                + std::imag( an[it] ) * std::imag( an[it] )
                                + std::real( bn[it] ) * std::real( bn[it] )
                                + std::imag( bn[it] ) * std::imag( bn[it] ) );
     }
     return Qsca * std::real(temp);
}


static double
C_Qext(iVec           an,
       iVec           bn,
       const double   x,
       const Vec      n)

{
     double Qsca = 2. / (x * x);
     complex128 temp;
     for(long unsigned int it = 0; it < an.size(); ++it)
     {
       temp += (2.*n[it]+1) * ( std::real( an[it] + an[it] ) );
     }
     return Qsca * std::real(temp);
}







static void
_S1S2(const double            index,
      const double            diameter,
      const double            wavelength,
      const double            nMedium,
      double                  *PhiPtr,
      const long unsigned int lenght,
      complex128*             s1s2)

{
    iVec an, bn;
    Vec n, n2;

    double m = index / nMedium,
           w = wavelength / nMedium,
           x = PI * diameter / w;



    int OrderMax = GetMaxOrder(x);

    std::tie(n, n2) = Arrange(1, OrderMax + 1);

    (x < 0.5) ? std::tie(an, bn) = LowFrequencyMie_ab(m, x) : std::tie(an, bn) = HighFrequencyMie_ab(m, x, OrderMax, n);

    iVec S1 = iVec(lenght),
         S2 = iVec(lenght),
         pin =iVec(OrderMax), taun = iVec(OrderMax);

    complex128 j (0., 1.0),
               *temp0 = &s1s2[0],
               *temp1 = &s1s2[lenght] ;

    for (long unsigned int i = 0; i < lenght; i++){

        std::tie(pin, taun) = MiePiTau(cos( PhiPtr[i]-PI/2 ), OrderMax);

        for (auto k = 0; k < OrderMax ; k++){
            *temp0 += n2[k] * ( an[k] * pin[k] +  bn[k] * taun[k] );
            *temp1 += n2[k] * ( an[k] * taun[k] + bn[k] * pin[k] ) ;

          }
    temp0++ ;
    temp1++ ;
    }


}







// -
