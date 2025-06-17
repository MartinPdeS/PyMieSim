#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include "../_bessel/bessel_subroutine.h"


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

      /**
       * @brief Computes the unstructured mode field.
       * @param x_coords X-coordinates of the points.
       * @param y_coords Y-coordinates of the points.
       * @return A vector of complex128 representing the unstructured mode field.
       */
      [[nodiscard]] std::vector<complex128> get_unstructured(std::vector<double> x_coords, std::vector<double> y_coords) const;


   private:
      /**
       * @brief Normalizes the coordinates to the range [-1, 1].
       * @param x_coords X-coordinates of the points.
       * @param y_coords Y-coordinates of the points.
       */
      void normalize_coordinates(std::vector<double> &x_coords, std::vector<double> &y_coords) const;

      /**
       * @brief Normalizes the field values to ensure they sum to 1.
       * @param field The field values to normalize.
       */
      void normalize_fields(std::vector<complex128> &field) const;

      /**
       * @brief Computes the Hermite-Gaussian mode field.
       * @param x_coords X-coordinates of the points.
       * @param y_coords Y-coordinates of the points.
       * @param wavelength Wavelength of the mode (default is 1.55).
       * @param waist_radius Waist radius of the mode (default is 0.3).
       * @param z Propagation distance (default is 0.0).
       * @return A vector of complex128 representing the Hermite-Gaussian mode field.
       */
      std::vector<complex128> get_LP_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords) const;

      /**
       * @brief Computes the Hermite-Gaussian mode field.
       * @param x_coords X-coordinates of the points.
       * @param y_coords Y-coordinates of the points.
       * @param wavelength Wavelength of the mode (default is 1.55).
       * @param waist_radius Waist radius of the mode (default is 0.3).
       * @param z Propagation distance (default is 0.0).
       * @return A vector of complex128 representing the Hermite-Gaussian mode field.
       */
      std::vector<complex128> get_HG_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords, double wavelength = 1.55, double waist_radius = 0.3, double z = 0.0) const ;

      /**
       * @brief Computes the Laguerre-Gaussian mode field.
       * @param x_coords X-coordinates of the points.
       * @param y_coords Y-coordinates of the points.
       * @param wavelength Wavelength of the mode (default is 1.55).
       * @param waist_radius Waist radius of the mode (default is 0.3).
       * @param z Propagation distance (default is 0.0).
       * @return A vector of complex128 representing the Laguerre-Gaussian mode field.
       */
      std::vector<complex128> get_LG_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords, double wavelength = 1.55, double waist_radius = 0.3, double z = 0.0) const ;

      /**
       * @brief Computes the non-circular unstructured mode field.
       * @param x_coords X-coordinates of the points.
       * @param y_coords Y-coordinates of the points.
       * @return A vector of complex128 representing the non-circular unstructured mode field.
       */
      std::vector<complex128> get_NC_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords) const;

      /**
       * @brief Computes the Hermite polynomial value for a given order and x.
       * @param n Order of the Hermite polynomial.
       * @param x Value at which to evaluate the polynomial.
       * @return The value of the Hermite polynomial at x.
       */
      double hermite_next(unsigned n, double x, double Hn, double Hnm1) const;

      /**
       * @brief Computes the Laguerre polynomial value for a given order and x.
       * @param n Order of the Laguerre polynomial.
       * @param x Value at which to evaluate the polynomial.
       * @return The value of the Laguerre polynomial at x.
       */
      double hermite_imp(unsigned n, double x) const;

      /**
       * @brief Computes the Laguerre polynomial value for a given order and x.
       * @param n Order of the Laguerre polynomial.
       * @param m Order of the associated Laguerre polynomial.
       * @param x Value at which to evaluate the polynomial.
       * @return The value of the Laguerre polynomial at x.
       */
      double laguerre_next(unsigned n, double x, double Ln, double Lnm1) const;

      /**
       * @brief Computes the Laguerre polynomial value for a given order and x.
       * @param n Order of the Laguerre polynomial.
       * @param m Order of the associated Laguerre polynomial.
       * @param x Value at which to evaluate the polynomial.
       * @return The value of the Laguerre polynomial at x.
       */
      double laguerre_next(unsigned n, unsigned l, double x, double Pl, double Plm1) const;

      /**
       * @brief Computes the Laguerre polynomial value for a given order and x.
       * @param n Order of the Laguerre polynomial.
       * @param x Value at which to evaluate the polynomial.
       * @return The value of the Laguerre polynomial at x.
       */
      double laguerre_imp(unsigned n, double x) const;

      /**
       * @brief Computes the Laguerre polynomial value for a given order and x.
       * @param n Order of the Laguerre polynomial.
       * @param m Order of the associated Laguerre polynomial.
       * @param x Value at which to evaluate the polynomial.
       * @return The value of the Laguerre polynomial at x.
       */
      double laguerre_imp(unsigned n, unsigned m, double x) const;


};
