#include "Functions.cpp"
#include <iostream>
namespace py = pybind11;

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;
#define j complex128(0.0,1.0)


 using namespace std;



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
S1S2(double   Index,
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
  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          SizeParam = PI * Diameter / Wavelength;

  int PhiLenght   = PhiBuffer.size,
      ThetaLenght = ThetaBuffer.size;

  iVec vec,
      __an = _an(SizeParam, Index, nMedium, MaxOrder),
      __bn = _bn(SizeParam, Index, nMedium, MaxOrder);


  Cndarray s1 = Cndarray(ThetaLenght),
           s2 = Cndarray(ThetaLenght);


  auto s1_data = s1.mutable_data(),
       s2_data = s2.mutable_data();

  for(auto l = 0; l < ThetaLenght; l++)
    {
      //tie(s1_data[l], s2_data[l]) = Expansion(BSC, __an, __bn, PhiPtr[l], ThetaPtr[l]) ;
      //s1_data[l] = 1.;
      //s2_data[l] = 1.;
    }

  return std::make_pair(s1, s2);
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

   double *PhiPtr   = (double *) PhiBuffer.ptr,
          *ThetaPtr = (double *) ThetaBuffer.ptr,
           SizeParam = PI * Diameter / Wavelength;

   int PhiLenght   = PhiBuffer.size,
       ThetaLenght = ThetaBuffer.size,
       Length      = BSCBuffer.shape[0],
       Width       = BSCBuffer.shape[1],
       n, m;

   double  *pin  = (double*) malloc(sizeof(double) * PhiLenght),
           *taun = (double*) malloc(sizeof(double) * PhiLenght),
           prefactor, n_f, m_f;

   complex128 *BSCPtr = (complex128 *) BSCBuffer.ptr,
              *an     = (complex128*) malloc(sizeof(complex128)*Length),
              *bn     = (complex128*) malloc(sizeof(complex128)*Length),
              TE, TM, _exp, __an, __bn, S1, S2;

   int first = (int)BSCPtr[0].real();

   CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

   Cndarray s1 = Cndarray(ThetaLenght*PhiLenght),
            s2 = Cndarray(ThetaLenght*PhiLenght);


   auto s1_data = s1.mutable_data(),
        s2_data = s2.mutable_data();

   for(auto p = 0; p < PhiLenght; p++)
     {
       //SymetricPinmTaunm(PhiPtr[p], pin, taun, Length);
       for(auto t = 0; t < ThetaLenght; t++)
       {
         for (auto b = 0; b < Length; b++)
         {
           S1 = 0.; S2 = 0.;
           n  = (int)BSCPtr[b].real();
           m  = (int)BSCPtr[b + Length].real();
           TE = BSCPtr[b + Length*2];
           TM = BSCPtr[b + Length*3];
           n_f = (double)n; m_f = (double)m;

           prefactor = (2.*n_f+1.)/( n_f* (n_f+1.) );
           _exp   = exp(j*m_f*ThetaPtr[b]);
           __an   = an[n-first+1];
           __bn   = bn[n-first+1];

           S1 += prefactor * ( ( j * __bn * TE * taun[b]) + (m_f * __an * TM * pin[b] ) ) * _exp;

           S2 += prefactor * ( ( j * m_f * __bn * TE * pin[b]) + (__an * TM * taun[b] ) ) * _exp;

           if (b == Length){s1_data[p] = S1; s2_data[p] = S2;}
         }
       }
     }

   return std::make_pair(s1, s2);
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

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
          SizeParam = PI * Diameter / Wavelength;

  int PhiLenght   = PhiBuffer.size,
      ThetaLenght = ThetaBuffer.size,
      Length      = BSCBuffer.shape[0],
      Width       = BSCBuffer.shape[1],
      n, m, w=0;

  double  *pin  = (double*) malloc(sizeof(double) * Length),
          *taun = (double*) malloc(sizeof(double) * Length),
          prefactor, n_f, m_f;

  complex128 *BSCPtr = (complex128 *) BSCBuffer.ptr,
             *an     = (complex128*) malloc(sizeof(complex128)*Length),
             *bn     = (complex128*) malloc(sizeof(complex128)*Length),
              TE, TM, _exp, __an, __bn, S1, S2;

  int first = (int)BSCPtr[0].real();

  CoefficientAnBn(SizeParam, Index, nMedium, MaxOrder, an, bn);

  Cndarray s1 = Cndarray(ThetaLenght*PhiLenght),
           s2 = Cndarray(ThetaLenght*PhiLenght);


  auto s1_data = s1.mutable_data(),
       s2_data = s2.mutable_data();

  for(auto p = 0; p < PhiLenght; p++)
    {
      SymetricPinmTaunm(Length, PhiPtr[p], pin, taun);
      for(auto t = 0; t < ThetaLenght; t++)
      {
        S1 = 0.; S2 = 0.;
        for (auto b = 0; b < 5; b++)
        {
          n  = (int)BSCPtr[b].real();
          m  = (int)BSCPtr[b + Length].real();
          TE = BSCPtr[b + Length*2];
          TM = BSCPtr[b + Length*3];
          n_f = (double)n; m_f = (double)m;

          prefactor = (2.*n_f+1.)/( n_f* (n_f+1.) );
          _exp   = exp( j * m_f * ThetaPtr[t] );
          __an   = an[n];
          __bn   = bn[n];


          S1 += prefactor * ( ( j * __bn * TE * taun[n]) + (m_f * __an * TM * pin[n] ) );// * _exp;

          cout << Length << "  " << n << "   TE=  " << PhiPtr[p] <<  "   TM=  "<< ThetaPtr[t] << m  << endl;

          S2 += prefactor * ( ( j * m_f * __bn * TE * pin[b]) + (__an * TM * taun[b] ) );// * _exp;

          //if (b == Length-1){s1_data[p+t] = S1; s2_data[p+t] = S2; cout<<S1<<endl;}
        }
        s1_data[w] = S1; s2_data[w] = S2;
        w++;
      }
    }

  return std::make_pair(s1, s2);
}



PYBIND11_MODULE(Sphere, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";



    module.def("S1S2",
               &S1S2,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("Polarization"),
               py::arg("E0"),
               py::arg("R"),
               py::arg("BSC"),
               py::arg("MaxOrder"),
               "Return S1");


     module.def("FieldsStructured",
                &FieldsStructured,
                py::arg("Index"),
                py::arg("Diameter"),
                py::arg("Wavelength"),
                py::arg("nMedium"),
                py::arg("Phi"),
                py::arg("Theta"),
                py::arg("Polarization"),
                py::arg("E0"),
                py::arg("R"),
                py::arg("BSC"),
                py::arg("MaxOrder"),
                "Return S1");


   module.def("an",
              &an,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return an");

   module.def("bn",
              &bn,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return bn");

   module.def("cn",
              &cn,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return cn");

   module.def("dn",
              &dn,
              py::arg("SizeParam"),
              py::arg("Index"),
              py::arg("nMedium"),
              py::arg("MaxOrder"),
              "Return dn");

}














//-
