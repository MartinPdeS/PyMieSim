double
GetQsca(complex128* an, complex128* bn, uint MaxOrder, double &SizeParam)
{
    double Qsca = 0.;

    for(uint it = 0; it < MaxOrder; ++it)
    {
         double n = (double) it + 1;
         Qsca += (2.* n + 1.) * ( pow( std::abs(an[it]), 2) + pow( std::abs(bn[it]), 2)  );

    }
    Qsca *= 2. / pow( SizeParam, 2.);
    return Qsca;
}


double
GetQext(complex128* an, complex128* bn, uint MaxOrder, double &SizeParam)
{
  double Qext = 0.;

  for(uint it = 0; it < MaxOrder; ++it)
  {
       double n = (double) it + 1;
       Qext += (2.* n + 1.) * std::real( an[it] + bn[it] );

  }
  Qext *= 2. / pow( SizeParam, 2.);
  return Qext;
}


double
Getg(complex128* an, complex128* bn, uint MaxOrder, double& SizeParam, double& Qsca)
{
    double g=0;

    for(uint it = 0; it < MaxOrder-1; ++it)
    {
      double n = (double) it + 1;
      g += ( n * (n + 2.) / (n + 1.) )            * std::real(an[it] * std::conj(an[it+1]) + bn[it] * std::conj(bn[it+1]) );
      g += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( an[it] * std::conj(bn[it]) );

    }

    g *= 4 / ( Qsca * pow(SizeParam, 2) );

    return g;
}


double
GetQback(complex128* an, complex128* bn, uint MaxOrder, double& SizeParam)
{
    complex128 Qback=0;

    for(uint it = 0; it < MaxOrder-1; ++it)
    {
      double n = (double) it + 1;
      Qback += (2. * n + 1) * pow(-1., n) * ( an[it] - bn[it] ) ;

    }
    Qback *= pow( std::abs(Qback), 2. ) / pow( SizeParam, 2. );
    return std::abs(Qback);
}

















// -
