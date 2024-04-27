#include <vector>
#include <cmath>
#include <definitions.cpp>
// #include <special_function.cpp>
#include <iostream>

// Factorial function using recursion
unsigned long long factorial(int n) {
    return (n == 0 || n == 1) ? 1 : n * factorial(n - 1);
}

// Function to calculate Bessel function of the first kind J_n(x) using series expansion
double bessel_J(int n, double x) {
    const int num_terms = 20;  // Number of terms in the series expansion
    double sum = 0.0;
    double x_half_pow_n = std::pow(x / 2.0, n);  // (x/2)^n part of the series

    for (int m = 0; m < num_terms; ++m) {
        double term = std::pow(-1.0, m) * std::pow(x / 2.0, 2 * m) / (factorial(m) * factorial(m + n));
        sum += term;
    }

    return x_half_pow_n * sum;
}

// Spherical Bessel function of the first kind j_n(x) using the relation to J_n(x)
double bessel_j(int n, double x) {
    if (x == 0) {
        return (n == 0) ? 1.0 : 0.0;  // j_0(0) = 1, j_n(0) = 0 for n > 0
    }
    double order = n + 0.5;
    double cylindrical_part = bessel_J(static_cast<int>(order), x);
    return std::sqrt(M_PI / (2 * x)) * cylindrical_part;
}

double bessel_jp(int n, double x)
{
    return 0.5 * (bessel_j(n - 1, x) - bessel_j(n + 1, x));
}

extern "C" {
    void jdzo(int& nt, int* n, int* m, double* zo);
}











void bjndd(int n, double x, std::vector<double>& bj, std::vector<double>& dj, std::vector<double>& fj) {
    if (x == 0) throw std::domain_error("x cannot be zero.");

    int m, mt, nt;
    double bs = 0.0, f0 = 0.0, f1 = 1e-35, f;

    // Determine the value of m based on the stability condition
    for (nt = 1; nt <= 900; ++nt) {
        mt = static_cast<int>(0.5 * log10(6.28 * nt) - nt * log10(1.36 * std::abs(x) / nt));
        if (mt > 20) break;
    }

    m = nt;
    // Recurrence for Bessel functions
    for (int k = m; k >= 0; --k) {
        f = 2.0 * (k + 1.0) * f1 / x - f0;
        if (k <= n) {
            bj[k] = f;
        }
        if (k % 2 == 0) {  // Check if k is even
            bs += 2.0 * f;
        }
        f0 = f1;
        f1 = f;
    }

    // Normalize Bessel functions
    for (int k = 0; k <= n; ++k) {
        bj[k] /= (bs - f);
    }

    // Calculate derivatives
    dj[0] = -bj[1];
    fj[0] = -bj[0] - dj[0] / x;
    for (int k = 1; k <= n; ++k) {
        dj[k] = bj[k-1] - k * bj[k] / x;
        fj[k] = (k * k / (x * x) - 1.0) * bj[k] - dj[k] / x;
    }
}

























std::vector<complex128> get_LP_mode_field(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int azimuthal_number,
   int radial_number)
{

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

   double u = 3.831705970207512315614435886308;

   #define N 40
   int nt = N; // Number of zeros
   int n[N], m[N];
   double zo[N];


   // Call the Fortran subroutine
   jdzo(nt, n, m, zo);

    // Printing all elements in the array
    for (int i = 0; i < N; i++)
        std::cout << "n: "<<n[i]<< "  zo[" << i << "] = " << zo[i] << std::endl;


   std::cout<<"U C++ " << u<<"\n";
   for (size_t i = 0; i < x_coords.size(); ++i)
   {
      double x = x_coords[i];
      double y = y_coords[i];

      double r = std::sqrt(x * x + y * y);
      double phi = std::atan2(y, x);

      // complex128 radial_part = bessel_j(azimuthal_number, u * r / core_radius);
      complex128 radial_part = bessel_J(azimuthal_number, r * u);

      double azimuthal_part = std::cos(azimuthal_number * phi);

      field[i] = radial_part * azimuthal_part;
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







// -