#pragma once

#include <vector>
#include <complex>
typedef std::complex<double> complex128;

namespace VSH
{
    namespace SPHERICAL
    {
        std::vector<complex128> compute_dn(double nmx, complex128 z); //Page 127 of BH

        std::tuple<std::vector<complex128>, std::vector<complex128>> MiePiTau(double mu, size_t max_order);

        void MiePiTau(double mu, size_t max_order, complex128 *pin, complex128 *taun);
    }

    namespace CYLINDRICAL
    {
      std::vector<complex128> compute_dn(double nmx, complex128 z); //Page 205 of BH
    }
}

