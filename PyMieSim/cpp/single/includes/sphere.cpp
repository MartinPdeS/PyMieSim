#pragma once

#include "single/headers/sphere.h"

typedef std::complex<double> complex128;

namespace SPHERE
{
    void Scatterer::compute_an_bn(){
        an.resize(max_order);
        bn.resize(max_order);

        complex128
            psi_n,
            chi_n,
            psi_1,
            chi_1,
            xi_n,
            xi_nm1,
            m = this->index / this->medium_index,
            mx = m * size_parameter,
            derivative_a, derivative_b;

        size_t nmx = std::max( max_order, (size_t) std::abs(mx) ) + 16;

        std::vector<complex128> Dn = VSH::SPHERICAL::compute_dn(nmx, mx);

        psi_1  = sin(size_parameter);
        chi_1 = cos(size_parameter);

        for (size_t order = 0; order < max_order; ++order)
        {
            // Calculate psi and chi (Riccati-Bessel functions)
            double nu = order + 1;
            psi_n = +size_parameter * compute_jn(nu, size_parameter);
            chi_n = -size_parameter * compute_yn(nu, size_parameter);

            // Complex Riccati-Bessel functions
            xi_n = psi_n - 1.0 * complex128(0, 1) * chi_n;
            xi_nm1 = psi_1 - 1.0 * complex128(0, 1) * chi_1;

            // Derivative of the Riccati-Bessel functions
            derivative_a = Dn[order + 1] / m + nu / size_parameter;
            derivative_b = Dn[order + 1] * m + nu / size_parameter;

            // Computation of the electric and magnetic multipole coefficients
            an[order] = (derivative_a * psi_n - psi_1) / (derivative_a * xi_n - xi_nm1);
            bn[order] = (derivative_b * psi_n - psi_1) / (derivative_b * xi_n - xi_nm1);

            psi_1 = psi_n;
            chi_1 = chi_n;
        }
    }

    void Scatterer::compute_cn_dn() {
        cn.resize(max_order);
        dn.resize(max_order);

        complex128
            x = size_parameter,
            m = this->index,
            z = m * x;

        size_t nmx = std::max( max_order, (size_t) std::abs(z) ) + 16;

        std::vector<complex128>
            Cnx = std::vector<complex128>(nmx),
            Cnn, jnx, jnmx, yx, hx, b1x, y1x, hn1x, ax, ahx, numerator,
            c_denominator, d_denominator;

        b1x.push_back( +sin(x) / x );
        y1x.push_back( -cos(x) / x );

        for (double i = nmx; i > 1; i--)
            Cnx[i-2] = i - z * z / (Cnx[i - 1] + i);

        for (size_t order = 0; order < max_order; order++)
        {
            Cnn.push_back(Cnx[order]);
            jnx.push_back(compute_jn(order + 1, x));

            jnmx.push_back(1.0 / (compute_jn(order + 1, z)));
            yx.push_back(compute_yn(order + 1, x));
            hx.push_back(jnx[order] + complex128(0, 1) * yx[order]);

            b1x.push_back(jnx[order]);
            y1x.push_back(yx[order]);
            hn1x.push_back(b1x[order] + complex128(0, 1) * y1x[order]);

            ax.push_back(x * b1x[order] - (complex128)(order + 1.0) * jnx[order]);
            ahx.push_back(x * hn1x[order] - (complex128)(order + 1.0) * hx[order]);

            numerator.push_back( jnx[order] * ahx[order] - hx[order] * ax[order] );
            c_denominator.push_back( ahx[order] - hx[order] * Cnn[order] );
            d_denominator.push_back( m * m * ahx[order] - hx[order] * Cnn[order] );
            cn[order] = jnmx[order] * numerator[order] / c_denominator[order] ;
            dn[order] = jnmx[order] * m * numerator[order] / d_denominator[order] ;
        }
    }

}

// -
