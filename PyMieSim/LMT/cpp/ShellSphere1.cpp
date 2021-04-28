#include <iostream>
#include <math.h>
#include <complex_bessel.h>

class SHELLSPHERE1: public BASE{

private:
  uint          MaxOrder;

  bool          Polarized;

  complex128    ShellIndex,
                CoreIndex;

  double        CoreDiameter,
                ShellDiameter,
                nMedium,
                xCore,
                xShell,
                Polarization,
                Wavelength,
                k,
                E0,
                Mu,
                MuScat;

    double&     Getk(){return this->k;};
    double&     GetPolarization(){return this->Polarization;};
    double&     GetE0(){return this->E0;};

    void        ComputeAnBn( complex128* an, complex128* bn, uint MaxOrder),
                LowFreqAnBn( complex128* an, complex128* bn),
                HighFreqAnBn(complex128* an, complex128* bn, uint MaxOrder),
                HighFreqCnDn(complex128* an, complex128* bn, uint MaxOrder);

    public:
      std::tuple<double, double, double, double, double, double, double> GetEfficiencies();

      Cndarray  An(uint MaxOrder),
                Bn(uint MaxOrder),
                Cn(uint MaxOrder),
                Dn(uint MaxOrder);


  SHELLSPHERE1(complex128 ShellIndex,
               complex128 CoreIndex,
               double     ShellDiameter,
               double     CoreDiameter,
               double     Wavelength,
               double     nMedium,
               double     Polarization,
               double     E0)
        {
          this->ShellDiameter = ShellDiameter;
          this->CoreDiameter  = CoreDiameter;
          this->ShellIndex    = ShellIndex;
          this->CoreIndex     = CoreIndex;
          this->nMedium       = nMedium;
          this->Wavelength    = Wavelength;
          this->E0            = E0;
          this->k             = 2 * PI / Wavelength;
          this->xShell        = k * ShellDiameter / 2;
          this->xCore         = k * CoreDiameter / 2;
          this->Polarization  = Polarization;
          this->MaxOrder      =  (uint) ( 2 + xShell + 4 * ( pow( xShell, 1./3. ) ) ) ;
          //this->SizeParam     = GetSizeParameter(Diameter, Wavelength, nMedium);
          this->Mu            = 1.0;
          this->MuScat        = 1.0;
        }

        ~SHELLSPHERE1(){  }

};

void
SHELLSPHERE1::ComputeAnBn(complex128* an, complex128* bn, uint MaxOrder)
{
  HighFreqAnBn(an, bn, MaxOrder);
}

Cndarray
SHELLSPHERE1::Bn(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return bn;
}


Cndarray
SHELLSPHERE1::An(uint MaxOrder)
{
  Cndarray     an         = Cndarray(MaxOrder),
               bn         = Cndarray(MaxOrder);

  complex128 * anPtr      = (complex128 *) an.request().ptr,
             * bnPtr      = (complex128 *) bn.request().ptr;

  this->HighFreqAnBn(anPtr, bnPtr, MaxOrder);
  return an;
}



void
SHELLSPHERE1::HighFreqAnBn(complex128* an, complex128* bn, uint MaxOrder)
{
  complex128 m = ShellIndex/CoreIndex,
             u = CoreIndex*xCore,
             v = ShellIndex*xCore,
             w = ShellIndex*xShell;

  complex128 sv = sqrt(0.5 * PI * v),
             sw = sqrt(0.5 * PI * w),
             sy = sqrt(0.5 * PI * xShell);

  uint mx   = (uint) std::max( abs( CoreIndex*xShell ), abs( ShellIndex*xShell ) );
  //int nmax = (int) ( 2. + xShell + 4. * ( pow(xShell, 1./3.) ) );
  int nmx  = (int) ( std::max( MaxOrder, mx ) + 16. )  ;

  iVec pv, pw, py, chv, chw, chy, p1y, ch1y, gsy, gs1y;

  p1y. push_back( sin(xShell) ) ;
  ch1y.push_back( cos(xShell) ) ;

  for (int i=0; i<MaxOrder+1; i++)
  {
    double nu = i + 1.5 ;
    pw.push_back( sw*F90Jn(nu,w) );
    pv.push_back( sv*F90Jn(nu,v) );
    py.push_back( sy*F90Jn(nu,xShell) );

    chv.push_back( -sv*F90Yn(nu,v) );
    chw.push_back( -sw*F90Yn(nu,w) );
    chy.push_back( -sy*F90Yn(nu,xShell) );

    p1y.push_back ( py[i]  );
    ch1y.push_back( chy[i] );
    gsy.push_back ( py[i]  - JJ * chy[i]  );
    gs1y.push_back ( p1y[i] - JJ * ch1y[i] );
  }

  iVec Du = iVec(nmx, 0.),
       Dv = iVec(nmx, 0.),
       Dw = iVec(nmx, 0.);

  for (int i = nmx-1; i > 1; i--)
  {
    Du[i-1] = (double)i / u -1. / (Du[i] + (double)i / u);
    Dv[i-1] = (double)i / v -1. / (Dv[i] + (double)i / v);
    Dw[i-1] = (double)i / w -1. / (Dw[i] + (double)i / w);
  }

  Du.erase(Du.begin());
  Dv.erase(Dv.begin());
  Dw.erase(Dw.begin());

  iVec uu, vv, fv, dns, gns, a1, b1;
  for (int i=0; i<MaxOrder; i++)
  {
    double n = (double) (i+1);
    uu.push_back ( m * Du[i] - Dv[i]  );
    vv.push_back ( Du[i] / m - Dv[i] );
    fv.push_back ( pv[i] / chv[i]    );
    dns.push_back( ( ( uu[i] * fv[i] / pw[i] ) / ( uu[i] * ( pw[i] - chw[i] * fv[i] ) + ( pw[i] / pv[i] ) / chv[i] ) ) + Dw[i] );
    gns.push_back( ( ( vv[i] * fv[i] / pw[i] ) / ( vv[i] * ( pw[i] - chw[i] * fv[i] ) + ( pw[i] / pv[i] ) / chv[i] ) ) + Dw[i] );
    a1.push_back ( dns[i] / ShellIndex + n / xShell );
    b1.push_back ( ShellIndex * gns[i] + n / xShell );
    an[i] = ( py[i] * a1[i] - p1y[i] ) / ( gsy[i] * a1[i] - gs1y[i] ) ;
    bn[i] = ( py[i] * b1[i] - p1y[i] ) / ( gsy[i] * b1[i] - gs1y[i] ) ;
  }



}

std::tuple<double, double, double, double, double, double, double>
SHELLSPHERE1::GetEfficiencies()
{
    double _nMedium = nMedium;

    complex128 _mCore   = CoreIndex  / _nMedium,
              _mShell  = ShellIndex / _nMedium;

    double _xCore   = PI * CoreDiameter  * _nMedium / Wavelength,
           _xShell  = PI * ShellDiameter * _nMedium / Wavelength;

    complex128 * an         = (complex128*) calloc(MaxOrder, sizeof(complex128)),
               * bn         = (complex128*) calloc(MaxOrder, sizeof(complex128));

    this->ComputeAnBn(an, bn, MaxOrder);

    double Qsca   = GetQsca(an, bn, MaxOrder, _xShell);
    double Qext   = GetQext(an, bn, MaxOrder, _xShell);
    std::cout<<Qext<<std::endl;
    double g      = Getg(an, bn, MaxOrder, _xShell, Qsca);
    double Qabs   = Qext - Qsca;
    double Qback  = GetQback(an, bn, MaxOrder, _xShell);
    double Qpr    = Qext - g * Qsca;
    double Qratio = Qback / Qsca;

    free(an);
    free(bn);
    return std::make_tuple(Qsca, Qext, Qabs, Qback, Qratio, g, Qpr);
}















// -
