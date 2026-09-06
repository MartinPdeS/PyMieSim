// amos_like_pure_cpp.hpp
#pragma once

#include <complex>
#include <cstdint>
#include <cmath>
#include <limits>

using complex128 = std::complex<double>;

#ifndef PI
#define PI (double)3.14159265358979323846264338
#endif

namespace amos_like_pure_cpp
{
    inline bool is_finite(complex128 z)
    {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    }

    inline complex128 complex_power(complex128 z, double a)
    {
        return std::exp(complex128(a, 0.0) * std::log(z));
    }

    // Lanczos approximation for complex log gamma
    inline complex128 complex_log_gamma(complex128 z)
    {
        static const double g = 7.0;
        static const double coefficients[] = {
            0.99999999999980993227684700473478,
            676.520368121885098567009190444019,
            -1259.13921672240287047156078755283,
            771.3234287776530788486528258894,
            -176.61502916214059906584551354,
            12.507343278686904814458936853,
            -0.13857109526572011689554707,
            9.984369578019570859563e-6,
            1.50563273514931155834e-7
        };

        if (z.real() < 0.5)
        {
            const complex128 pi(PI, 0.0);
            return std::log(pi) - std::log(std::sin(pi * z)) - complex_log_gamma(complex128(1.0, 0.0) - z);
        }

        z -= complex128(1.0, 0.0);

        complex128 series_sum(coefficients[0], 0.0);
        for (int i = 1; i < 9; ++i)
        {
            series_sum += complex128(coefficients[i], 0.0) / (z + complex128(static_cast<double>(i), 0.0));
        }

        const complex128 t = z + complex128(g + 0.5, 0.0);
        const complex128 half_log_two_pi(0.91893853320467274178032973640562, 0.0); // 0.5*log(2*pi)

        return half_log_two_pi + (z + complex128(0.5, 0.0)) * std::log(t) - t + std::log(series_sum);
    }

    inline complex128 complex_gamma(complex128 z)
    {
        return std::exp(complex_log_gamma(z));
    }

    inline void set_outputs_from_complex(complex128 value, double& out_real, double& out_imag, int& nz, int& ierr)
    {
        nz = 0;

        if (!is_finite(value))
        {
            out_real = std::numeric_limits<double>::quiet_NaN();
            out_imag = std::numeric_limits<double>::quiet_NaN();
            ierr = 3;
            return;
        }

        out_real = value.real();
        out_imag = value.imag();
        ierr = 0;
    }

    inline complex128 bessel_j_series(double order, complex128 z, int& ierr)
    {
        const double epsilon = 1e-15;
        const int maximum_iterations = 20000;

        const complex128 half_z = z / complex128(2.0, 0.0);
        complex128 term = complex_power(half_z, order) / complex_gamma(complex128(order + 1.0, 0.0));

        complex128 series_sum = term;

        for (int k = 1; k < maximum_iterations; ++k)
        {
            const double k_as_double = static_cast<double>(k);
            term *= complex128(-1.0, 0.0) * half_z * half_z;
            term /= complex128(k_as_double, 0.0);
            term /= complex128(order + k_as_double, 0.0);

            series_sum += term;

            if (std::abs(term) < epsilon * std::abs(series_sum))
            {
                ierr = 0;
                return series_sum;
            }
        }

        ierr = 2;
        return series_sum;
    }

    inline complex128 bessel_i_series(double order, complex128 z, int& ierr)
    {
        const double epsilon = 1e-15;
        const int maximum_iterations = 20000;

        const complex128 half_z = z / complex128(2.0, 0.0);
        complex128 term = complex_power(half_z, order) / complex_gamma(complex128(order + 1.0, 0.0));

        complex128 series_sum = term;

        for (int k = 1; k < maximum_iterations; ++k)
        {
            const double k_as_double = static_cast<double>(k);
            term *= half_z * half_z;
            term /= complex128(k_as_double, 0.0);
            term /= complex128(order + k_as_double, 0.0);

            series_sum += term;

            if (std::abs(term) < epsilon * std::abs(series_sum))
            {
                ierr = 0;
                return series_sum;
            }
        }

        ierr = 2;
        return series_sum;
    }

    inline complex128 bessel_y_from_j(double order, complex128 z, int& ierr)
    {
        const complex128 pi(PI, 0.0);
        const complex128 nu(order, 0.0);

        const complex128 sin_pi_nu = std::sin(pi * nu);
        const complex128 cos_pi_nu = std::cos(pi * nu);

        if (std::abs(sin_pi_nu) < 1e-12)
        {
            // Near integer order: this formula is numerically unstable.
            // We use a tiny perturbation. This is not AMOS accurate but avoids division by near zero.
            const double small_perturbation = 1e-10;
            return bessel_y_from_j(order + small_perturbation, z, ierr);
        }

        int ierr_j_nu = 0;
        int ierr_j_minus_nu = 0;

        const complex128 j_nu = bessel_j_series(order, z, ierr_j_nu);
        const complex128 j_minus_nu = bessel_j_series(-order, z, ierr_j_minus_nu);

        ierr = (ierr_j_nu != 0) ? ierr_j_nu : ierr_j_minus_nu;

        // Y_nu(z) = ( J_nu(z) cos(pi nu) - J_{-nu}(z) ) / sin(pi nu)
        return (j_nu * cos_pi_nu - j_minus_nu) / sin_pi_nu;
    }

    inline complex128 bessel_k_from_i(double order, complex128 z, int& ierr)
    {
        const complex128 pi(PI, 0.0);
        const complex128 nu(order, 0.0);

        const complex128 sin_pi_nu = std::sin(pi * nu);

        if (std::abs(sin_pi_nu) < 1e-12)
        {
            const double small_perturbation = 1e-10;
            return bessel_k_from_i(order + small_perturbation, z, ierr);
        }

        int ierr_i_nu = 0;
        int ierr_i_minus_nu = 0;

        const complex128 i_nu = bessel_i_series(order, z, ierr_i_nu);
        const complex128 i_minus_nu = bessel_i_series(-order, z, ierr_i_minus_nu);

        ierr = (ierr_i_nu != 0) ? ierr_i_nu : ierr_i_minus_nu;

        // K_nu(z) = (pi/2) * (I_{-nu}(z) - I_nu(z)) / sin(pi nu)
        return (pi / complex128(2.0, 0.0)) * (i_minus_nu - i_nu) / sin_pi_nu;
    }
}

// Bessel function of the first kind.
inline void zbesj_wrap(
    double zr,
    double zi,
    double order,
    int32_t kode,
    int32_t N,
    double* cyr,
    double* cyi,
    int32_t* nz,
    int32_t* ierr)
{
    int nz_local = 0;
    int ierr_local = 0;

    if (kode < 1 || kode > 2 || N != 1 || !std::isfinite(order))
    {
        *cyr = 0.0;
        *cyi = 0.0;
        *nz = 0;
        *ierr = 1;
        return;
    }

    const complex128 z(zr, zi);
    complex128 value = amos_like_pure_cpp::bessel_j_series(order, z, ierr_local);

    if (kode == 2)
    {
        value *= std::exp(-std::abs(z.imag()));
    }

    double out_r = 0.0;
    double out_i = 0.0;
    amos_like_pure_cpp::set_outputs_from_complex(value, out_r, out_i, nz_local, ierr_local);

    *cyr = out_r;
    *cyi = out_i;
    *nz = static_cast<int32_t>(nz_local);
    *ierr = static_cast<int32_t>(ierr_local);
}

// Bessel function of the second kind.
inline void zbesy_wrap(
    double zr,
    double zi,
    double order,
    int32_t kode,
    int32_t N,
    double* cyr,
    double* cyi,
    int32_t* nz,
    double* /*cwrkr*/,
    double* /*cwrki*/,
    int32_t* ierr)
{
    int nz_local = 0;
    int ierr_local = 0;

    if (kode < 1 || kode > 2 || N != 1 || !std::isfinite(order))
    {
        *cyr = 0.0;
        *cyi = 0.0;
        *nz = 0;
        *ierr = 1;
        return;
    }

    const complex128 z(zr, zi);

    if (zr == 0.0 && zi == 0.0)
    {
        *cyr = 0.0;
        *cyi = 0.0;
        *nz = 0;
        *ierr = 1;
        return;
    }

    complex128 value = amos_like_pure_cpp::bessel_y_from_j(order, z, ierr_local);

    if (kode == 2)
    {
        value *= std::exp(-std::abs(z.imag()));
    }

    double out_r = 0.0;
    double out_i = 0.0;
    amos_like_pure_cpp::set_outputs_from_complex(value, out_r, out_i, nz_local, ierr_local);

    *cyr = out_r;
    *cyi = out_i;
    *nz = static_cast<int32_t>(nz_local);
    *ierr = static_cast<int32_t>(ierr_local);
}

// Modified Bessel function of the first kind.
inline void zbesi_wrap(double zr, double zi, double order, int32_t kode, int32_t N,
                       double& cyr, double& cyi, int32_t& nz, int32_t& ierr)
{
    int nz_local = 0;
    int ierr_local = 0;

    if (kode < 1 || kode > 2 || N != 1 || !std::isfinite(order))
    {
        cyr = 0.0;
        cyi = 0.0;
        nz = 0;
        ierr = 1;
        return;
    }

    const complex128 z(zr, zi);
    complex128 value = amos_like_pure_cpp::bessel_i_series(order, z, ierr_local);

    if (kode == 2)
    {
        value *= std::exp(-std::abs(z.real()));
    }

    double out_r = 0.0;
    double out_i = 0.0;
    amos_like_pure_cpp::set_outputs_from_complex(value, out_r, out_i, nz_local, ierr_local);

    cyr = out_r;
    cyi = out_i;
    nz = static_cast<int32_t>(nz_local);
    ierr = static_cast<int32_t>(ierr_local);
}

// Modified Bessel function of the second kind.
inline void zbesk_wrap(double zr, double zi, double order, int32_t kode, int32_t N,
                       double& cyr, double& cyi, int32_t& nz, int32_t& ierr)
{
    int nz_local = 0;
    int ierr_local = 0;

    if (kode < 1 || kode > 2 || N != 1 || !std::isfinite(order))
    {
        cyr = 0.0;
        cyi = 0.0;
        nz = 0;
        ierr = 1;
        return;
    }

    const complex128 z(zr, zi);

    if (zr == 0.0 && zi == 0.0)
    {
        cyr = 0.0;
        cyi = 0.0;
        nz = 0;
        ierr = 1;
        return;
    }

    complex128 value = amos_like_pure_cpp::bessel_k_from_i(order, z, ierr_local);

    if (kode == 2)
    {
        value *= std::exp(z);
    }

    double out_r = 0.0;
    double out_i = 0.0;
    amos_like_pure_cpp::set_outputs_from_complex(value, out_r, out_i, nz_local, ierr_local);

    cyr = out_r;
    cyi = out_i;
    nz = static_cast<int32_t>(nz_local);
    ierr = static_cast<int32_t>(ierr_local);
}

// Hankel functions of both kinds.
inline void zbesh_wrap(
    double zr,
    double zi,
    double order,
    int32_t kode,
    int32_t kind,
    int32_t N,
    double* cyr,
    double* cyi,
    int32_t* nz,
    int32_t* ierr)
{
    if (!cyr || !cyi || !nz || !ierr)
    {
        return;
    }

    int nz_local = 0;
    int ierr_local = 0;

    *nz = 0;

    if (N != 1 || (kind != 1 && kind != 2) || (kode != 1 && kode != 2) || !std::isfinite(order))
    {
        *cyr = 0.0;
        *cyi = 0.0;
        *ierr = 1;
        return;
    }

    const complex128 z(zr, zi);

    int ierr_j = 0;
    int ierr_y = 0;

    const complex128 J = amos_like_pure_cpp::bessel_j_series(order, z, ierr_j);
    const complex128 Y = amos_like_pure_cpp::bessel_y_from_j(order, z, ierr_y);

    ierr_local = (ierr_j != 0) ? ierr_j : ierr_y;

    const complex128 i_unit(0.0, 1.0);
    complex128 value = (kind == 1) ? (J + i_unit * Y) : (J - i_unit * Y);

    if (kode == 2)
    {
        const double factor = static_cast<double>(3 - 2 * kind); // kind=1 -> +1, kind=2 -> -1
        value *= std::exp(-factor * i_unit * z);
    }

    double out_r = 0.0;
    double out_i = 0.0;
    amos_like_pure_cpp::set_outputs_from_complex(value, out_r, out_i, nz_local, ierr_local);

    *cyr = out_r;
    *cyi = out_i;
    *nz  = static_cast<int32_t>(nz_local);
    *ierr = static_cast<int32_t>(ierr_local);
}