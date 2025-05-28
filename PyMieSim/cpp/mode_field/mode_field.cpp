#include "mode_field.h"


#define PI (double)3.14159265358979323846264338
typedef std::complex<double> complex128;


// ________________________________ COMMON ___________________________________
// ________________________________ STD ___________________________________

void ModeField::normalize_coordinates(std::vector<double> &x_coords, std::vector<double> &y_coords) const {
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
}

void ModeField::normalize_fields(std::vector<complex128> &field) const {
   double norm = 0.0;
   for (const complex128& f : field)
      norm += std::norm(f);

   norm = std::sqrt(norm);

   if (norm != 0) {  // Avoid division by zero
      for (complex128& f : field)
         f /= norm;
   }
}


// ________________________________ LP ___________________________________

[[nodiscard]] std::vector<complex128>
ModeField::get_LP_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords) const {

   size_t azimuthal_number = mode_id.number_0;
   size_t radial_number = mode_id.number_1;

   std::vector<complex128> field(x_coords.size());

   // Normalize the coordinates
   this->normalize_coordinates(x_coords, y_coords);

   std::vector<double> rj0(radial_number), rj1(radial_number), ry0(radial_number), ry1(radial_number);

   // Calculate Bessel zeros for the radial and azimuthal components
   bessel_zeros(azimuthal_number, radial_number, &rj0[0], &rj1[0], &ry0[0], &ry1[0]);

   double x, y, r, phi, azimuthal_part;

   for (size_t i = 0; i < x_coords.size(); ++i) {
      x = x_coords[i];
      y = y_coords[i];

      r = std::sqrt(x * x + y * y);
      phi = std::atan2(y, x);

      complex128 radial_part = bessel_J(azimuthal_number, r * rj0[radial_number-1]);

      azimuthal_part = std::cos(azimuthal_number * phi);

      field[i] = radial_part * azimuthal_part;
   }

   // Normalization to L2 norm of 1
   this->normalize_fields(field);

   return field;
}



// ________________________________ HG ___________________________________

double ModeField::hermite_next(unsigned n, double x, double Hn, double Hnm1) const {
   return (2 * x * Hn - 2 * n * Hnm1);
}


double ModeField::hermite_imp(unsigned n, double x) const {
   double p0 = 1;
   double p1 = 2 * x;

   if(n == 0)
      return p0;

   unsigned c = 1;

   while(c < n) {
      std::swap(p0, p1);
      p1 = this->hermite_next(c, x, p0, p1);
      ++c;
   }
   return p1;
}


// Helper function to calculate the Hermite-Gaussian mode field amplitude
std::vector<complex128> ModeField::get_HG_unstructured(std::vector<double>& x_coords, std::vector<double>& y_coords, double wavelength, double waist_radius, double z) const {

   size_t x_number = this->mode_id.number_0;
   size_t y_number = this->mode_id.number_1;

   double k = 2 * PI / wavelength;  // Wave number
   double w0 = waist_radius;  // Beam waist

   std::vector<complex128> field(x_coords.size());

   // Normalize the coordinates
   this->normalize_coordinates(x_coords, y_coords);

   // Calculate the beam width at distance z
   double w = w0 * std::sqrt(1 + (z * wavelength / (PI * w0 * w0)) * (z * wavelength / (PI * w0 * w0)));

   // Radius of curvature of the beam's wavefront at z
   double R = (z == 0) ? std::numeric_limits<double>::infinity() :
   z * (1 + (PI * w0 * w0 / (z * wavelength)) * (PI * w0 * w0 / (z * wavelength)));

   // Gouy phase shift at z
   double gouy_phase = std::atan(z * PI / (wavelength * w0 * w0));

   // Process each coordinate
   for (size_t i = 0; i < x_coords.size(); ++i) {
      double x = x_coords[i];
      double y = y_coords[i];

      // Hermite polynomial factors
      double Hn = this->hermite_imp(x_number, std::sqrt(2) * x / w);
      double Hm = this->hermite_imp(y_number, std::sqrt(2) * y / w);

      // Amplitude calculation
      double amplitude = Hn * Hm * std::exp(-((x * x + y * y) / (w * w)));

      // Phase calculation including Gouy phase and spherical phase factor
      double phase = -k * ((x * x + y * y) / (2 * R)) + (x_number + y_number) * gouy_phase;

      // Combine amplitude and phase into a complex field
      field[i] = complex128(amplitude * std::cos(phase), amplitude * std::sin(phase));
   }

   // Normalization to L2 norm of 1
   this->normalize_fields(field);

   return field;
}





// ________________________________ LG ___________________________________


double ModeField::laguerre_next(unsigned n, double x, double Ln, double Lnm1) const
{
   return ((2 * n + 1 - x) * Ln - n * Lnm1) / (n + 1);
}


double ModeField::laguerre_next(unsigned n, unsigned l, double x, double Pl, double Plm1) const
{
   return ((2 * n + l + 1 - x) * Pl - (n + l) * Plm1) / (n + 1);
}


double ModeField::laguerre_imp(unsigned n, double x) const
{
   double p0 = 1;
   double p1 = 1 - x;

   if(n == 0)
      return p0;

   unsigned c = 1;

   while(c < n) {
      std::swap(p0, p1);
      p1 = laguerre_next(c, x, p0, p1);
      ++c;
   }
   return p1;
}


double ModeField::laguerre_imp(unsigned n, unsigned m, double x) const
{
   // Special cases:
   if(m == 0)
      return laguerre_imp(n, x);

   double p0 = 1;

   if(n == 0)
      return p0;

   double p1 = m + 1 - x;

   unsigned c = 1;

   while(c < n) {
      std::swap(p0, p1);
      p1 = laguerre_next(c, m, x, p0, p1);
      ++c;
   }
   return p1;
}


std::vector<complex128> ModeField::get_LG_unstructured(std::vector<double> &x_coords, std::vector<double> &y_coords, double wavelength, double waist_radius, double z) const {

   size_t azimuthal_number = mode_id.number_0;
   size_t radial_number = mode_id.number_1;


   double k = 2 * PI / wavelength;  // Wave number
   double w0 = waist_radius;  // Beam waist

   std::vector<complex128> field(x_coords.size());

   // Normalize the coordinates
   this->normalize_coordinates(x_coords, y_coords);

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
      double L_pl = laguerre_imp(azimuthal_number, radial_number, 2 * r * r / (w * w));
      double amplitude = std::pow(std::sqrt(2) * r / w, radial_number) * L_pl * std::exp(-r * r / (w * w));

      // Phase factor
      double phase = radial_number * theta - k * r * r / (2 * R) + (2 * azimuthal_number + radial_number + 1) * gouy_phase;

      field[i] = amplitude * std::cos(phase); // store the real part of the field
    }

   // Normalization to L2 norm of 1
   this->normalize_fields(field);


   return field;
}


// ________________________________ NC ___________________________________


std::vector<complex128> ModeField::get_NC_unstructured(std::vector<double> &x_coords, std::vector<double> &) const {
   std::vector<complex128> output(x_coords.size(), 1.0);

   return output;
}

