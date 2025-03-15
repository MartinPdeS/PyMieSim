#include "utils/numpy_interface.cpp"
#include <vector>
#include <cmath>
#include <complex>

#define PI (double)3.14159265358979323846264338
typedef std::complex<double> complex128;


double laguerre_next(unsigned n, double x, double Ln, double Lnm1)
{
   return ((2 * n + 1 - x) * Ln - n * Lnm1) / (n + 1);
}


double laguerre_next(unsigned n, unsigned l, double x, double Pl, double Plm1)
{
   return ((2 * n + l + 1 - x) * Pl - (n + l) * Plm1) / (n + 1);
}


double laguerre_imp(unsigned n, double x)
{
   double p0 = 1;
   double p1 = 1 - x;

   if(n == 0)
      return p0;

   unsigned c = 1;

   while(c < n)
   {
      std::swap(p0, p1);
      p1 = laguerre_next(c, x, p0, p1);
      ++c;
   }
   return p1;
}


double laguerre_imp(unsigned n, unsigned m, double x)
{
   // Special cases:
   if(m == 0)
      return laguerre_imp(n, x);

   double p0 = 1;

   if(n == 0)
      return p0;

   double p1 = m + 1 - x;

   unsigned c = 1;

   while(c < n)
   {
      std::swap(p0, p1);
      p1 = laguerre_next(c, m, x, p0, p1);
      ++c;
   }
   return p1;
}


std::vector<complex128> get_LG_mode_field(
    std::vector<double> x_coords,
    std::vector<double> y_coords,
    int azimuthal_number,
    int radial_number,
    double wavelength = 1.55,
    double waist_radius = 0.3,
    double z = 0.0)
{
   double k = 2 * PI / wavelength;  // Wave number
   double w0 = waist_radius;  // Beam waist

   std::vector<complex128> field(x_coords.size());

   // Normalize the coordinates
   double max_norm = 0.0;
   for (size_t i = 0; i < x_coords.size(); ++i) {
      double norm = std::sqrt(x_coords[i] * x_coords[i] + y_coords[i] * y_coords[i]);
      if (norm > max_norm)
         max_norm = norm;

   }

   if (max_norm != 0) {  // Avoid division by zero
      for (size_t i = 0; i < x_coords.size(); ++i) {
         x_coords[i] /= max_norm;
         y_coords[i] /= max_norm;
      }
   }

   // Convert to polar coordinates and calculate field values
   for (size_t i = 0; i < x_coords.size(); ++i) {
        double x = x_coords[i];
        double y = y_coords[i];
        double r = std::sqrt(x * x + y * y);
        double theta = std::atan2(y, x);

      // Beam parameters at z
      double w = w0 * std::sqrt(1 + (z * wavelength / (PI * w0 * w0)) * (z * wavelength / (PI * w0 * w0)));
      double R = (z == 0) ? std::numeric_limits<double>::infinity() : z * (1 + (PI * w0 * w0 / (z * wavelength)) * (PI * w0 * w0 / (z * wavelength)));
      double gouy_phase = std::atan(z * PI / (wavelength * w0 * w0));

      // Laguerre polynomial
      double L_pl = laguerre_imp(azimuthal_number, std::abs(radial_number), 2 * r * r / (w * w));
      double amplitude = std::pow(std::sqrt(2) * r / w, std::abs(radial_number)) * L_pl * std::exp(-r * r / (w * w));

      // Phase factor
      double phase = radial_number * theta - k * r * r / (2 * R) + (2 * azimuthal_number + std::abs(radial_number) + 1) * gouy_phase;

      field[i] = amplitude * std::cos(phase); // store the real part of the field
    }

   // Normalization to L2 norm of 1
   double norm = 0.0;
   for (const complex128& f : field)
       norm += std::norm(f);

   norm = std::sqrt(norm);

   for (complex128& f : field)
      f /= norm;


   return field;
}

pybind11::array_t<complex128> get_LG_mode_field_py(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int azimuthal_number,
   int radial_number)
{
   return _vector_to_numpy(
      get_LG_mode_field(x_coords, y_coords, azimuthal_number, radial_number),
      {x_coords.size()}
   );
}




