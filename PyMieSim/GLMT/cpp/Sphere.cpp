#include <iostream>
#define j complex128(0.0,1.0)


namespace Sphere
{

  iVec
  _an(double SizeParam, double Index, double nMedium, int MaxOrder)
  {
    double alpha = SizeParam,
           beta = alpha * Index,
           MuSp = 1.,
           Mu = 1.,
           M = Index/nMedium;

    complex128 numerator, denominator;
    iVec _an;
    for (auto order = 1; order < MaxOrder+1; order++)
    {
      numerator = MuSp * Psi(order, alpha) * Psi_p(order, beta)  - Mu * M * Psi_p(order, alpha) * Psi(order, beta);
      denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
      _an.push_back(numerator/denominator);
    }
    return _an;
  }


  iVec
  _bn(double SizeParam, double Index, double nMedium, int MaxOrder)
  {
    double alpha = SizeParam,
           beta = alpha * Index,
           MuSp = 1.,
           Mu = 1.,
           M = Index/nMedium;

    complex128 numerator, denominator;
    iVec _bn;
    for (auto order = 1; order < MaxOrder+1; order++)
    {
      numerator = Mu * M * Psi(order, alpha) * Psi_p(order, beta) - MuSp * Psi_p(order, alpha) * Psi(order, beta);
      denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp  * Xi_p(order, alpha) * Psi(order, beta);
      _bn.push_back(numerator/denominator);
    }
    return _bn;
  }



  iVec
  _cn(double SizeParam, double Index, double nMedium, int MaxOrder)
  {
    double alpha = SizeParam,
           beta = alpha * Index,
           MuSp = 1.,
           Mu = 1.,
           M = Index/nMedium;

    complex128 numerator, denominator;
    iVec _cn;
    for (auto order = 1; order < MaxOrder+1; order++)
    {
      numerator = M * MuSp * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
      denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
      _cn.push_back(numerator/denominator);
    }
    return _cn;
  }


  iVec
  _dn(double SizeParam, double Index, double nMedium, int MaxOrder)
  {

    double alpha = SizeParam,
           beta = alpha * Index,
           MuSp = 1.,
           Mu = 1.,
           M = Index/nMedium;

    complex128 numerator, denominator;
    iVec _dn;
    for (auto order = 1; order < MaxOrder+1; order++)
    {
      numerator = Mu * M*M * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
      denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp * Xi_p(order, alpha) * Psi(order, beta);
      _dn.push_back(numerator/denominator);
    }
    return _dn;
  }




void
CoefficientAnBn(const double &SizeParam,
           const double &Index,
           const double &nMedium,
           const int    &MaxOrder,
           complex128   *an,
           complex128   *bn)
{
  double alpha = SizeParam,
         beta  = alpha * Index, 
         MuSp  = 1.,
         Mu    = 1.,
         M     = Index/nMedium;

  complex128 numerator, denominator, PsiAlpha, PsiBeta, PsiPBeta, PsiPAlpha, XiAlpha, XiPAlpha;

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    PsiAlpha  = Psi(order, alpha);
    PsiBeta   = Psi(order, beta);
    PsiPBeta  = Psi_p(order, beta);
    PsiPAlpha = Psi_p(order, alpha);
    XiAlpha   = Xi(order, alpha);
    XiPAlpha  = Xi_p(order, alpha);

    numerator = MuSp * PsiAlpha * PsiPBeta  - Mu * M * PsiPAlpha * PsiBeta;
    denominator = MuSp * XiAlpha * PsiPBeta - Mu * M * XiPAlpha * PsiBeta;
    an[order] = numerator/denominator;

    numerator = Mu * M * PsiAlpha * PsiPBeta - MuSp * PsiPAlpha * PsiBeta;
    denominator = Mu * M * XiAlpha * PsiPBeta - MuSp  * XiPAlpha * PsiBeta;
    bn[order] = numerator/denominator;
  }
}



Cndarray
an(double SizeParam, double Index, double nMedium, int MaxOrder)
{
  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  Cndarray _an = Cndarray(MaxOrder);
  auto an_data = _an.mutable_data();

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = MuSp * Psi(order, alpha) * Psi_p(order, beta)  - Mu * M * Psi_p(order, alpha) * Psi(order, beta);
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    an_data[order-1] = numerator/denominator;
  }
  return _an;
}


Cndarray
bn(double SizeParam, double Index, double nMedium, int MaxOrder)
{

  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  Cndarray _bn = Cndarray(MaxOrder);
  auto bn_data = _bn.mutable_data();

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M * Psi(order, alpha) * Psi_p(order, beta) - MuSp * Psi_p(order, alpha) * Psi(order, beta);
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp  * Xi_p(order, alpha) * Psi(order, beta);
    bn_data[order-1] = numerator/denominator;
  }
  return _bn;
}



Cndarray
cn(double SizeParam, double Index, double nMedium, int MaxOrder)
{

  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  Cndarray _cn = Cndarray(MaxOrder);
  auto cn_data = _cn.mutable_data();

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = M * MuSp * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = MuSp * Xi(order, alpha) * Psi_p(order, beta) - Mu * M * Xi_p(order, alpha) * Psi(order, beta);
    cn_data[order-1] = numerator/denominator;
  }
  return _cn;
}


Cndarray
dn(double SizeParam, double Index, double nMedium, int MaxOrder)
{

  double alpha = SizeParam,
         beta = alpha * Index,
         MuSp = 1.,
         Mu = 1.,
         M = Index/nMedium;

  complex128 numerator, denominator;
  Cndarray _dn = Cndarray(MaxOrder);
  auto dn_data = _dn.mutable_data();

  for (auto order = 1; order < MaxOrder+1; order++)
  {
    numerator = Mu * M*M * ( Xi(order, alpha) * Psi_p(order, alpha) - Xi_p(order, alpha) * Psi(order, alpha) );
    denominator = Mu * M * Xi(order, alpha) * Psi_p(order, beta) - MuSp * Xi_p(order, alpha) * Psi(order, beta);
    dn_data[order-1] = numerator/denominator;
  }
  return _dn;
}





std::pair<Cndarray, Cndarray>
FieldsUnstructured(double   Index,
                   double   Diameter,
                   double   Wavelength,
                   double   nMedium,
                   ndarray  Phi,
                   ndarray  Theta,
                   double   Polarization,
                   double   E0,
                   double   R,
                   Cndarray BSC,
                   int      MaxOrder)
 {
   py::buffer_info PhiBuffer     = Phi.request();
   py::buffer_info ThetaBuffer   = Theta.request();
   py::buffer_info BSCBuffer = BSC.request();

   int ThetaLenght = ThetaBuffer.size,
       Length      = BSCBuffer.shape[0],
       n           = 0,
       m           = 0;

   double *PhiPtr    = (double *) PhiBuffer.ptr,
          *ThetaPtr  = (double *) ThetaBuffer.ptr,
           SizeParam = PI * Diameter / Wavelength,
           k         = 2. * PI / Wavelength,
          *pin       = (double*) malloc(sizeof(complex128)*Length),
          *taun      = (double*) malloc(sizeof(complex128)*Length),
           prefactor = 0.;

   complex128 *an     = (complex128*) malloc(sizeof(complex128)*Length),
              *bn     = (complex128*) malloc(sizeof(complex128)*Length),
              *BSCPtr = (complex128 *) BSCBuffer.ptr,
               _exp   = 0.,
               S1     = 0.,
               S2     = 0.,
               TE     = 0.,
               TM     = 0.,
               propagator = E0 / (k * R) * exp(-j*k*R);

   CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

   Cndarray s1 = Cndarray(ThetaLenght),
            s2 = Cndarray(ThetaLenght);

   auto s1_data = s1.mutable_data(),
        s2_data = s2.mutable_data();

   for(auto l = 0; l < ThetaLenght; l++)
     {
       SymetricPinmTaunm(Length, PhiPtr[l]-PI/2, pin, taun);
       S1 = 0.; S2=0.;
       for (auto b = 0; b < Length; b++)
         {
             n  = (int)BSCPtr[b].real();
             m  = (int)BSCPtr[b + Length].real();
             TE = BSCPtr[b + Length*2];
             TM = BSCPtr[b + Length*3];
             prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );
             _exp   = exp( j * (double)m * ThetaPtr[l] + Polarization );

             S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

             S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;

        }
        s1_data[l] = S1 * propagator;
        s2_data[l] = S2 * propagator;

     }

    free(an);
    free(bn);
    free(pin);
    free(taun);
   return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
 }



 std::pair<Cndarray, Cndarray>
 FieldsStructured(double   Index,
                  double   Diameter,
                  double   Wavelength,
                  double   nMedium,
                  ndarray  Phi,
                  ndarray  Theta,
                  double   Polarization,
                  double   E0,
                  double   R,
                  Cndarray BSC,
                  int      MaxOrder)
  {
    py::buffer_info PhiBuffer     = Phi.request();
    py::buffer_info ThetaBuffer   = Theta.request();
    py::buffer_info BSCBuffer = BSC.request();

    int PhiLenght   = PhiBuffer.size,
        ThetaLenght = ThetaBuffer.size,
        Length      = BSCBuffer.shape[0],
        n           = 0,
        m           = 0,
        index       = 0;

    double *PhiPtr    = (double *) PhiBuffer.ptr,
           *ThetaPtr  = (double *) ThetaBuffer.ptr,
            k         = 2. * PI / Wavelength,
            SizeParam = PI * Diameter / Wavelength,
           *pin       = (double*)    malloc(sizeof(complex128)*Length),
           *taun      = (double*)   malloc(sizeof(complex128)*Length),
            prefactor = 0.;

    complex128 *an         = (complex128*) malloc(sizeof(complex128)*Length),
               *bn         = (complex128*) malloc(sizeof(complex128)*Length),
               *BSCPtr     = (complex128 *) BSCBuffer.ptr,
                _exp       = 0.,
                S1         = 0.,
                S2         = 0.,
                TE         = 0.,
                TM         = 0.,
                propagator = E0 / (k * R) * exp(-j*k*R);

    CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

    Cndarray s1 = Cndarray(ThetaLenght*PhiLenght),
             s2 = Cndarray(ThetaLenght*PhiLenght);



    complex128 *s1_data = s1.mutable_data(),
               *s2_data = s2.mutable_data();

    for(auto p = 0; p < PhiLenght; p++)
      {
        SymetricPinmTaunm(Length, PhiPtr[p]-PI/2, pin, taun);

        for(auto t = 0; t < ThetaLenght; t++)
          {
              S1 = 0.; S2=0.;
              for (auto b = 0; b < Length; b++)
                {
                    n  = (int)BSCPtr[b].real();
                    m  = (int)BSCPtr[b + Length].real();
                    TE = BSCPtr[b + Length*2];
                    TM = BSCPtr[b + Length*3];
                    prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );
                    _exp   = exp( j * (double)m * ThetaPtr[t] + Polarization );

                    S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

                    S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;


               }

               s1_data[index] = S1 * propagator;
               s2_data[index] = S2 * propagator;
               index++;
          }

      }

     free(an);
     free(bn);
     free(pin);
     free(taun);

     s1.resize({PhiLenght,ThetaLenght});
     s2.resize({PhiLenght,ThetaLenght});
     return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
  }





std::pair<Cndarray, Cndarray>
FieldsUnstructuredUnpolarized(double   Index,
                               double   Diameter,
                               double   Wavelength,
                               double   nMedium,
                               ndarray  Phi,
                               ndarray  Theta,
                               double   E0,
                               double   R,
                               Cndarray BSC,
                               int      MaxOrder)
 {
   py::buffer_info PhiBuffer     = Phi.request();
   py::buffer_info ThetaBuffer   = Theta.request();
   py::buffer_info BSCBuffer = BSC.request();

   int ThetaLenght = ThetaBuffer.size,
       Length      = BSCBuffer.shape[0],
       n           = 0,
       m           = 0;

   double *PhiPtr    = (double *) PhiBuffer.ptr,
           SizeParam = PI * Diameter / Wavelength,
           k         = 2. * PI / Wavelength,
          *pin       = (double*) malloc(sizeof(complex128)*Length),
          *taun      = (double*) malloc(sizeof(complex128)*Length),
           prefactor = 0.;

   complex128 *an     = (complex128*) malloc(sizeof(complex128)*Length),
              *bn     = (complex128*) malloc(sizeof(complex128)*Length),
              *BSCPtr = (complex128 *) BSCBuffer.ptr,
               _exp   = 1./sqrt(2.),
               S1     = 0.,
               S2     = 0.,
               TE     = 0.,
               TM     = 0.,
               propagator = E0 / (k * R) * exp(-j*k*R);

   CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

   Cndarray s1 = Cndarray(ThetaLenght),
            s2 = Cndarray(ThetaLenght);

   auto s1_data = s1.mutable_data(),
        s2_data = s2.mutable_data();

   for(auto l = 0; l < ThetaLenght; l++)
     {
       SymetricPinmTaunm(Length, PhiPtr[l]-PI/2, pin, taun);
       S1 = 0.; S2=0.;
       for (auto b = 0; b < Length; b++)
         {
             n  = (int)BSCPtr[b].real();
             m  = (int)BSCPtr[b + Length].real();
             TE = BSCPtr[b + Length*2];
             TM = BSCPtr[b + Length*3];
             prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );

             S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

             S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;

        }
        s1_data[l] = S1 * propagator;
        s2_data[l] = S2 * propagator;

     }

    free(an);
    free(bn);
    free(pin);
    free(taun);
   return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
 }



 std::pair<Cndarray, Cndarray>
 FieldsStructuredUnpolarized(double   Index,
                              double   Diameter,
                              double   Wavelength,
                              double   nMedium,
                              ndarray  Phi,
                              ndarray  Theta,
                              double   E0,
                              double   R,
                              Cndarray BSC,
                              int      MaxOrder)
  {
    py::buffer_info PhiBuffer     = Phi.request();
    py::buffer_info ThetaBuffer   = Theta.request();
    py::buffer_info BSCBuffer = BSC.request();

    int PhiLenght   = PhiBuffer.size,
        ThetaLenght = PhiBuffer.size,
        Length      = BSCBuffer.shape[0],
        n           = 0,
        m           = 0,
        index       = 0;

    double *PhiPtr    = (double *) PhiBuffer.ptr,
            k         = 2. * PI / Wavelength,
            SizeParam = PI * Diameter / Wavelength,
           *pin       = (double*)    malloc(sizeof(complex128)*Length),
           *taun      = (double*)   malloc(sizeof(complex128)*Length),
            prefactor = 0.;

    complex128 *an         = (complex128*) malloc(sizeof(complex128)*Length),
               *bn         = (complex128*) malloc(sizeof(complex128)*Length),
               *BSCPtr     = (complex128 *) BSCBuffer.ptr,
                _exp       = 1./sqrt(2.),
                S1         = 0.,
                S2         = 0.,
                TE         = 0.,
                TM         = 0.,
                propagator = E0 / (k * R) * exp(-j*k*R);

    CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

    Cndarray s1 = Cndarray(ThetaLenght*PhiLenght),
             s2 = Cndarray(ThetaLenght*PhiLenght);



    auto s1_data = s1.mutable_data(),
         s2_data = s2.mutable_data();

    for(auto p = 0; p < PhiLenght; p++)
      {
        SymetricPinmTaunm(Length, PhiPtr[p]-PI/2, pin, taun);

        for(auto t = 0; t < ThetaLenght; t++)
          {
              S1 = 0.; S2=0.;
              for (auto b = 0; b < Length; b++)
                {
                    n  = (int)BSCPtr[b].real();
                    m  = (int)BSCPtr[b + Length].real();
                    TE = BSCPtr[b + Length*2];
                    TM = BSCPtr[b + Length*3];
                    prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );

                    S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

                    S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;

               }
               s1_data[index] = S1 * propagator;
               s2_data[index] = S2 * propagator;
               index++;
          }

      }

     free(an);
     free(bn);
     free(pin);
     free(taun);
     s1.resize({PhiLenght,ThetaLenght});
     s2.resize({PhiLenght,ThetaLenght});
     return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
  }






std::pair<Cndarray, Cndarray>
S1S2Unstructured(double   Index,
                 double   Diameter,
                 double   Wavelength,
                 double   nMedium,
                 ndarray  Phi,
                 ndarray  Theta,
                 double   Polarization,
                 Cndarray BSC,
                 int      MaxOrder)
  {
 py::buffer_info PhiBuffer     = Phi.request();
 py::buffer_info ThetaBuffer   = Theta.request();
 py::buffer_info BSCBuffer = BSC.request();

 int ThetaLenght = ThetaBuffer.size,
     Length      = BSCBuffer.shape[0],
     n           = 0,
     m           = 0;

 double *PhiPtr    = (double *) PhiBuffer.ptr,
        *ThetaPtr  = (double *) ThetaBuffer.ptr,
         SizeParam = PI * Diameter / Wavelength,
        *pin       = (double*)    malloc(sizeof(complex128)*Length),
        *taun      = (double*)   malloc(sizeof(complex128)*Length),
         prefactor = 0.;

 complex128 *an     = (complex128*) malloc(sizeof(complex128)*Length),
            *bn     = (complex128*) malloc(sizeof(complex128)*Length),
            *BSCPtr = (complex128 *) BSCBuffer.ptr,
             _exp   = 0.,
             S1     = 0.,
             S2     = 0.,
             TE     = 0.,
             TM     = 0.;

 CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

 Cndarray s1 = Cndarray(ThetaLenght),
          s2 = Cndarray(ThetaLenght);

 auto s1_data = s1.mutable_data(),
      s2_data = s2.mutable_data();

 for(auto l = 0; l < ThetaLenght; l++)
   {
     SymetricPinmTaunm(Length, PhiPtr[l]-PI/2, pin, taun);
     S1 = 0.; S2=0.;
     for (auto b = 0; b < Length; b++)
       {
           n  = (int)BSCPtr[b].real();
           m  = (int)BSCPtr[b + Length].real();
           TE = BSCPtr[b + Length*2];
           TM = BSCPtr[b + Length*3];
           prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );
           _exp   = exp( j * (double)m * ThetaPtr[l] + Polarization );

           S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

           S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;

      }
      s1_data[l] = S1;
      s2_data[l] = S2;

   }

  free(an);
  free(bn);
  free(pin);
  free(taun);
 return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
}



std::pair<Cndarray, Cndarray>
S1S2Structured(double   Index,
               double   Diameter,
               double   Wavelength,
               double   nMedium,
               ndarray  Phi,
               ndarray  Theta,
               double   Polarization,
               Cndarray BSC,
               int      MaxOrder)
{
  py::buffer_info PhiBuffer     = Phi.request();
  py::buffer_info ThetaBuffer   = Theta.request();
  py::buffer_info BSCBuffer = BSC.request();

  int PhiLenght   = PhiBuffer.size,
      ThetaLenght = ThetaBuffer.size,
      Length      = BSCBuffer.shape[0],
      n           = 0,
      m           = 0,
      index       = 0;

  double *PhiPtr    = (double *) PhiBuffer.ptr,
         *ThetaPtr  = (double *) ThetaBuffer.ptr,
          SizeParam = PI * Diameter / Wavelength,
         *pin       = (double*)    malloc(sizeof(complex128)*Length),
         *taun      = (double*)   malloc(sizeof(complex128)*Length),
          prefactor = 0.;

  complex128 *an     = (complex128*) malloc(sizeof(complex128)*Length),
             *bn     = (complex128*) malloc(sizeof(complex128)*Length),
             *BSCPtr = (complex128 *) BSCBuffer.ptr,
              _exp   = 0.,
              S1     = 0.,
              S2     = 0.,
              TE     = 0.,
              TM     = 0.;

  CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

  Cndarray s1 = Cndarray(ThetaLenght*PhiLenght),
           s2 = Cndarray(ThetaLenght*PhiLenght);



  auto s1_data = s1.mutable_data(),
       s2_data = s2.mutable_data();

  for(auto p = 0; p < PhiLenght; p++)
    {
      SymetricPinmTaunm(Length, PhiPtr[p]-PI/2, pin, taun);

      for(auto t = 0; t < ThetaLenght; t++)
        {
            S1 = 0.; S2=0.;
            for (auto b = 0; b < Length; b++)
              {
                  n  = (int)BSCPtr[b].real();
                  m  = (int)BSCPtr[b + Length].real();
                  TE = BSCPtr[b + Length*2];
                  TM = BSCPtr[b + Length*3];
                  prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );
                  _exp   = exp( j * (double)m * ThetaPtr[t] + Polarization );

                  S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

                  S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;

             }
             s1_data[index] = S1;
             s2_data[index] = S2;
             index++;
        }

    }

   free(an);
   free(bn);
   free(pin);
   free(taun);
   s1.resize({PhiLenght,ThetaLenght});
   s2.resize({PhiLenght,ThetaLenght});

   return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
}








std::pair<Cndarray, Cndarray>
S1S2UnstructuredUnpolarized(double  Index,
                           double   Diameter,
                           double   Wavelength,
                           double   nMedium,
                           ndarray  Phi,
                           ndarray  Theta,
                           Cndarray BSC,
                           int      MaxOrder)
  {
 py::buffer_info PhiBuffer     = Phi.request();
 py::buffer_info ThetaBuffer   = Theta.request();
 py::buffer_info BSCBuffer = BSC.request();

 int ThetaLenght = ThetaBuffer.size,
     Length      = BSCBuffer.shape[0],
     n           = 0,
     m           = 0;

 double *PhiPtr    = (double *) PhiBuffer.ptr,
         SizeParam = PI * Diameter / Wavelength,
        *pin       = (double*)    malloc(sizeof(complex128)*Length),
        *taun      = (double*)   malloc(sizeof(complex128)*Length),
         prefactor = 0.;

 complex128 *an     = (complex128*) malloc(sizeof(complex128)*Length),
            *bn     = (complex128*) malloc(sizeof(complex128)*Length),
            *BSCPtr = (complex128 *) BSCBuffer.ptr,
             _exp   = 1./sqrt(2.),
             S1     = 0.,
             S2     = 0.,
             TE     = 0.,
             TM     = 0.;

 CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

 Cndarray s1 = Cndarray(ThetaLenght),
          s2 = Cndarray(ThetaLenght);

 auto s1_data = s1.mutable_data(),
      s2_data = s2.mutable_data();

 for(auto l = 0; l < ThetaLenght; l++)
   {
     SymetricPinmTaunm(Length, PhiPtr[l]-PI/2, pin, taun);
     S1 = 0.; S2=0.;
     for (auto b = 0; b < Length; b++)
       {
           n  = (int)BSCPtr[b].real();
           m  = (int)BSCPtr[b + Length].real();
           TE = BSCPtr[b + Length*2];
           TM = BSCPtr[b + Length*3];
           prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );

           S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

           S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;

      }
      s1_data[l] = S1;
      s2_data[l] = S2;

   }

  free(an);
  free(bn);
  free(pin);
  free(taun);
 return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
}



std::pair<Cndarray, Cndarray>
S1S2StructuredUnpolarized(double   Index,
                         double   Diameter,
                         double   Wavelength,
                         double   nMedium,
                         ndarray  Phi,
                         ndarray  Theta,
                         Cndarray BSC,
                         int      MaxOrder)
{
  py::buffer_info PhiBuffer     = Phi.request();
  py::buffer_info ThetaBuffer   = Theta.request();
  py::buffer_info BSCBuffer     = BSC.request();

  int PhiLenght   = PhiBuffer.size,
      ThetaLenght = ThetaBuffer.size,
      Length      = BSCBuffer.shape[0],
      n           = 0,
      m           = 0,
      index       = 0;

  double *PhiPtr    = (double *) PhiBuffer.ptr,
          SizeParam = PI * Diameter / Wavelength,
         *pin       = (double*)    malloc(sizeof(complex128)*Length),
         *taun      = (double*)   malloc(sizeof(complex128)*Length),
          prefactor = 0.;

  complex128 *an     = (complex128*) malloc(sizeof(complex128)*Length),
             *bn     = (complex128*) malloc(sizeof(complex128)*Length),
             *BSCPtr = (complex128 *) BSCBuffer.ptr,
              _exp   = 1./sqrt(2.),
              S1     = 0.,
              S2     = 0.,
              TE     = 0.,
              TM     = 0.;

  CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

  Cndarray s1 = Cndarray(ThetaLenght*PhiLenght),
           s2 = Cndarray(ThetaLenght*PhiLenght);



  auto s1_data = s1.mutable_data(),
       s2_data = s2.mutable_data();

  for(auto p = 0; p < PhiLenght; p++)
    {
      SymetricPinmTaunm(Length, PhiPtr[p]-PI/2, pin, taun);

      for(auto t = 0; t < ThetaLenght; t++)
        {
            S1 = 0.; S2=0.;
            for (auto b = 0; b < Length; b++)
              {
                  n         = (int)BSCPtr[b].real();
                  m         = (int)BSCPtr[b + Length].real();
                  TE        = BSCPtr[b + Length*2];
                  TM        = BSCPtr[b + Length*3];
                  prefactor = (2. * (double)n + 1.) / ( (double)n * ( (double)n + 1. ) );

                  S1 += prefactor*(j * bn[n] * TE * taun[n] + (double)m * an[n] * TM * pin[n]) * _exp;

                  S2 += prefactor*(j * (double)m * bn[n] * TE * pin[n] +  an[n] * TM * taun[n]) * _exp;

             }
             s1_data[index] = S1;
             s2_data[index] = S2;
             index++;
        }

    }

   free(an);
   free(bn);
   free(pin);
   free(taun);
   s1.resize({PhiLenght,ThetaLenght});
   s2.resize({PhiLenght,ThetaLenght});
   return std::make_pair(s1.attr("transpose")(), s2.attr("transpose")());
}




}








//-
