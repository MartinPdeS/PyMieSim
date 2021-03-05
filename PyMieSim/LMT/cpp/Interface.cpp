#include "Functions.cpp"

namespace py = pybind11;

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef py::array_t<double> ndarray;
typedef py::array_t<complex128> Cndarray;

#define j complex128(0.0,1.0)



std::tuple<double, double, double>
Efficiencies(const double  m, const double  x)

{
    iVec an, bn;

    double Qsca, Qext, Qabs;

    Vec n, n2;

    int OrderMax = GetMaxOrder(x);

    std::tie(n, n2) = Arrange(1, OrderMax + 1);

    (x < 0.5) ? std::tie(an, bn) = LowFrequencyMie_ab(m, x) : std::tie(an, bn) = HighFrequencyMie_ab(m, x, OrderMax, n);

    Qsca = C_Qsca(an, bn, x, n);

    Qext = C_Qext(an, bn, x, n);

    Qabs = Qext - Qsca;

    return std::make_tuple(Qsca, Qext, Qabs);
}


std::tuple<Cndarray, Cndarray>
S1S2(const double            index,
     const double            diameter,
     const double            wavelength,
     const double            nMedium,
     ndarray                 Phi)

{
    py::buffer_info PhiBuffer   = Phi.request();

    double *PhiPtr   = (double *) PhiBuffer.ptr,
           m = index / nMedium,
           w = wavelength / nMedium,
           x = PI * diameter / w;

    int lenght = PhiBuffer.size,
        OrderMax = GetMaxOrder(x);

    iVec pin = iVec(OrderMax),
         taun = iVec(OrderMax),
         an,
         bn;

    Vec n, n2;

    complex128 temp0=0., temp1=0.;

    Cndarray S1 = ndarray(lenght,0),
             S2 = ndarray(lenght,0);

    auto S1_data = S1.mutable_data(),
         S2_data = S2.mutable_data();

    std::tie(n, n2) = Arrange(1, OrderMax + 1);

    (x < 0.5) ? std::tie(an, bn) = LowFrequencyMie_ab(m, x) : std::tie(an, bn) = HighFrequencyMie_ab(m, x, OrderMax, n);

    for (auto i = 0; i < lenght; i++){

        std::tie(pin, taun) = MiePiTau(cos( PhiPtr[i]-PI/2. ), OrderMax);

        for (auto k = 0; k < OrderMax ; k++){
            temp0 += n2[k] * ( an[k] * pin[k] +  bn[k] * taun[k] );
            temp1 += n2[k] * ( an[k] * taun[k] + bn[k] * pin[k]  );
          }
        S1_data[i] = temp0;
        S2_data[i] = temp1;

        temp0 = 0.; temp1=0.;

    }

    return std::make_tuple(S1,S2);
}


std::pair<Cndarray, Cndarray>
Fields(double     index,
       double     diameter,
       double     wavelength,
       double     nMedium,
       ndarray    Phi,
       ndarray    Theta,
       double     Polarization,
       double     E0,
       double     R,
       int        Lenght)
{

  py::buffer_info PhiBuffer   = Phi.request();
  py::buffer_info ThetaBuffer = Theta.request();

  py::array_t<complex128> result0 = py::array_t<complex128>(Phi.size());
  py::array_t<complex128> result1 = py::array_t<complex128>(Phi.size());

  py::buffer_info buf0 = result0.request();
  py::buffer_info buf1 = result1.request();

  double *PhiPtr   = (double *) PhiBuffer.ptr,
         *ThetaPtr = (double *) ThetaBuffer.ptr,
         k = 2. * PI / wavelength,
         temp0 ;

  complex128 *ptr0 = (complex128 *) buf0.ptr,
             *ptr1 = (complex128 *) buf1.ptr,
             temp2,
             propagator = E0 / (k * R) * exp(-j*k*R),
             *s1s2 = (complex128*) calloc(2 * Lenght , sizeof(complex128));

  ndarray Parallel, Perpendicular;

  _S1S2(index,
       diameter,
       wavelength,
       nMedium,
       PhiPtr,
       Lenght,
       s1s2);

  for (auto k=0; k < Lenght; k++ )
  {
    temp0 = ThetaPtr[k] ;
    ptr0[k] = j* propagator * s1s2[k] * (complex128) abs(cos(temp0 + Polarization));
    ptr1[k] =- propagator * s1s2[k + Lenght] * (complex128) abs(sin(temp0 + Polarization));

  }

  free(s1s2);

  return std::make_pair(result0,result1);
}


PYBIND11_MODULE(Sphere, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) c++ binding module for light scattering from a spherical scatterer";

    module.def("Fields",
               &Fields,
               py::arg("Index"),
               py::arg("Diameter"),
               py::arg("Wavelength"),
               py::arg("nMedium"),
               py::arg("Phi"),
               py::arg("Theta"),
               py::arg("Polarization"),
               py::arg("E0"),
               py::arg("R"),
               py::arg("Lenght"),
               "Compute the scattering far-field for a spherical scatterer");


     module.def("S1S2",
                &S1S2,
                py::arg("Index"),
                py::arg("Diameter"),
                py::arg("Wavelength"),
                py::arg("nMedium"),
                py::arg("Phi"),
                "Compute the scattering coefficient S1 & S2");


      module.def("Efficiencies",
                 &Efficiencies,
                 py::arg("Index"),
                 py::arg("SizeParameter"),
                 "Compute the scattering efficiencies");

}







// -
