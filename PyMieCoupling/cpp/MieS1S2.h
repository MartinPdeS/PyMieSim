
#include <iostream>
#include <complex>
#include <vector>
#include <tuple>

typedef std::vector<std::complex<double>> iVec;
typedef std::complex<double> complex128;


void Sum1(iVec &);


std::tuple<std::vector<complex128> , std::vector<complex128>> LowFrequencyMie_ab(double m,double x);
