#pragma once

typedef std::complex<double> complex128;

namespace VSH{
    namespace SPHERICAL {

        std::vector<complex128> compute_dn(double nmx, complex128 z) //Page 127 of BH
        {
          std::vector<complex128> Dn(nmx, 0.0);

          for (double n = nmx - 1.; n > 1.; n--)
              Dn[n-1] = n/z - ( 1. / (Dn[n] + n/z) );

           return Dn;
        }

        inline std::tuple<std::vector<complex128>, std::vector<complex128>> MiePiTau(double mu, size_t max_order)
        {
          std::vector<complex128> pin, taun;
          pin.reserve(max_order);
          taun.reserve(max_order);

          pin.push_back( 1. );
          pin.push_back( 3. * mu );

          taun.push_back( mu );
          taun.push_back( 3.0 * cos(2. * acos(mu) ) );

          for (size_t order = 2; order < max_order; order++)
              {
               pin.push_back( ( (2. * (double)order + 1.) * mu * pin[order - 1] - ((double)order + 1.) * pin[order - 2] ) / (double)order );

               taun.push_back( ((double)order + 1.) * mu * pin[order] - ((double)order + 2.) * pin[order - 1] );
             }

          return std::make_tuple(pin, taun);
        }

        inline void MiePiTau(double mu, size_t max_order, complex128 *pin, complex128 *taun)
        {
          pin[0] = 1.;
          pin[1] = 3. * mu;

          taun[0] = mu;
          taun[1] = 3.0 * cos(2. * acos(mu) );

          for (size_t order = 2; order < max_order; order++)
              {
               pin[order] = ( (2. * (double)order + 1.) * mu * pin[order - 1] - ((double)order + 1.) * pin[order - 2] ) / (double)order;

               taun[order] = ((double)order + 1.) * mu * pin[order] - ((double)order + 2.) * pin[order - 1];
             }
        }
    }

    namespace CYLINDRICAL {

      std::vector<complex128> compute_dn(double nmx, complex128 z) //Page 205 of BH
      {
        std::vector<complex128> Dn(nmx, 0.0);

        for (double n = nmx - 1; n > 0; n--)
            Dn[n-1] = n/z - ( 1. / (Dn[n] + n/z) );

         return Dn;
      }

    }
}

