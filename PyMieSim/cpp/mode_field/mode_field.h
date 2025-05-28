#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include <pybind11/numpy.h>
#include "../_bessel/bessel_subroutine.h"
#include "utils/numpy_interface.h"


#define PI (double)3.14159265358979323846264338
typedef std::complex<double> complex128;


struct ModeID {
   std::string mode_family;
   int number_0;
   int number_1;
};

class ModeField {
   ModeID mode_id;

   public:
      ModeField() = default;
      ModeField(const ModeID &mode_id) : mode_id(mode_id) {}

      [[nodiscard]] std::vector<complex128> get_unstructured(std::vector<double> x_coords, std::vector<double> y_coords) const {
         if (this->mode_id.mode_family == "LP")
            return this->get_LP_unstructured(x_coords, y_coords);
         if (this->mode_id.mode_family == "HG")
            return this->get_HG_unstructured(x_coords, y_coords);
         if (this->mode_id.mode_family == "LG")
            return this->get_LG_unstructured(x_coords, y_coords);
         if (this->mode_id.mode_family == "NC")
            return this->get_NC_unstructured(x_coords, y_coords);

         throw std::runtime_error("Invalid mode family");
      };


   private:
      void normalize_coordinates(std::vector<double> &x_coords, std::vector<double> &y_coords) const;
      void normalize_fields(std::vector<complex128> &field) const;

      std::vector<complex128> get_LP_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords) const;
      std::vector<complex128> get_HG_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords, double wavelength = 1.55, double waist_radius = 0.3, double z = 0.0) const ;
      std::vector<complex128> get_LG_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords, double wavelength = 1.55, double waist_radius = 0.3, double z = 0.0) const ;
      std::vector<complex128> get_NC_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords) const;

      double hermite_next(unsigned n, double x, double Hn, double Hnm1) const;
      double hermite_imp(unsigned n, double x) const;

      double laguerre_next(unsigned n, double x, double Ln, double Lnm1) const;
      double laguerre_next(unsigned n, unsigned l, double x, double Pl, double Plm1) const;
      double laguerre_imp(unsigned n, double x) const;
      double laguerre_imp(unsigned n, unsigned m, double x) const;


};
