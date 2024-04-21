#pragma once

#include "definitions.cpp"
#include "fibonacci_mesh.cpp"


namespace VSH{
    namespace SPHERICAL {

        std::vector<complex128> compute_dn(double &&nmx, complex128 &z) //Page 127 of BH
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

          double n = 0;
          for (unsigned int i = 2; i < max_order; i++)
              {
               n = (double)i;

               pin.push_back( ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n );

               taun.push_back( (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1] );
             }

          return std::make_tuple(pin, taun);
        }

        inline void MiePiTau(double mu, unsigned int max_order, complex128 *pin, complex128 *taun)
        {
          pin[0] = 1.;
          pin[1] = 3. * mu;

          taun[0] = mu;
          taun[1] = 3.0 * cos(2. * acos(mu) );

          double n = 0;
          for (unsigned int i = 2; i < max_order; i++)
              {
               n = (double)i;

               pin[i] = ( (2. * n + 1.) * mu * pin[i-1] - (n + 1.) * pin[i-2] ) / n;

               taun[i] = (n + 1.) * mu * pin[i] - (n + 2.) * pin[i-1];
             }
        }




        SphericalCoordinate RDir(double &theta, double &phi)
        {
            SphericalCoordinate Output;
            Output.r     = sin(theta)*cos(phi);
            Output.phi   = sin(theta)*sin(phi),
            Output.theta = cos(theta);
            return Output;
        }

        SphericalCoordinate ThetaDir(double &theta, double &phi)
        {
            SphericalCoordinate Output;
            Output.r     = cos(theta)*cos(phi);
            Output.phi   = cos(theta)*sin(phi),
            Output.theta = -sin(theta);
            return Output;
        }

        SphericalCoordinate PhiDir(double &theta, double &phi)
        {
            SphericalCoordinate Output;
            Output.r     = -sin(phi);
            Output.phi   = cos(phi),
            Output.theta = 0;
            return Output;
        }
    }

    namespace CYLINDRICAL {

      std::vector<complex128> compute_dn(double &&nmx, complex128 &z) //Page 205 of BH
      {
        std::vector<complex128> Dn(nmx, 0.0);

        for (double n = nmx - 1; n > 0; n--)
            Dn[n-1] = n/z - ( 1. / (Dn[n] + n/z) );

         return Dn;
      }

    }
}

