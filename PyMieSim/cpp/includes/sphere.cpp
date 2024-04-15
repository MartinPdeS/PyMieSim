#pragma once

#include "sphere.h"

namespace SPHERE
{
    using complex128 = std::complex<double>;

    void Scatterer::compute_max_order(size_t max_order) {
        if (max_order == 0) {
            this->max_order = static_cast<size_t>(2 + this->size_parameter + 4 * std::cbrt(this->size_parameter)) + 16;
        } else {
            this->max_order = max_order;
        }
    }

    void Scatterer::compute_size_parameter() {
        this->size_parameter = PI * this->diameter / this->source.wavelength;
    }

    void Scatterer::compute_area() {
        this->area = PI * std::pow(this->diameter / 2.0, 2);
    }

    void Scatterer::compute_an_bn(){
        an.resize(max_order);
        bn.resize(max_order);

        complex128
            x = size_parameter,
            m = this->index,
            mx = m * x,
            _da, _db, _gsx, _gs1x, _px, _chx, _p1x, _ch1x, _p2x, _ch2x;

        size_t nmx = std::max( max_order, (size_t) std::abs(mx) ) + 16;

        std::vector<complex128> Dn = VSH::SPHERICAL::compute_dn(nmx, mx);

        double n;

        _p1x  = sin(x);
        _ch1x = cos(x);

        for (size_t i = 1; i < max_order + 1; ++i)
        {
            n = (double) i;
            _px =  x * compute_jn(n, x);
            _chx = -x * compute_yn(n, x);

            _p2x = _px;
            _ch2x = _chx;

            _gsx =  _px  - 1.0 * JJ * _chx;
            _gs1x =  _p1x - 1.0 * JJ * _ch1x;

            _da = Dn[i] / m + n / x;
            _db = Dn[i] * m + n / x;

            an[i-1] = (_da * _px - _p1x) / (_da * _gsx - _gs1x);
            bn[i-1] = (_db * _px - _p1x) / (_db * _gsx - _gs1x);

            _p1x  = _p2x;
            _ch1x = _ch2x;
        }
    }

    void Scatterer::compute_cn_dn(){
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
        {
            Cnx[i-2] = i - z*z/(Cnx[i-1] + i);
        }

        for (size_t i = 0; i < max_order; i++)
        {
            double n = (double) i;

            Cnn.push_back(Cnx[i]);
            jnx.push_back(compute_jn( n+1, x ));

            jnmx.push_back(1. / ( compute_jn(n+1, z )));
            yx.push_back(compute_yn(n+1, x ));
            hx.push_back(jnx[i] + JJ * yx[i]);

            b1x.push_back(jnx[i]);
            y1x.push_back(yx[i]);
            hn1x.push_back(b1x[i] + JJ * y1x[i]);

            ax.push_back(x * b1x[i] - ( n+1 ) * jnx[i]);
            ahx.push_back(x * hn1x[i] - ( n+1 ) * hx[i]);

            numerator.push_back( jnx[i] * ahx[i] - hx[i] * ax[i] );
            c_denominator.push_back( ahx[i] - hx[i] * Cnn[i] );
            d_denominator.push_back( m * m * ahx[i] - hx[i] * Cnn[i] );
            cn[i] = jnmx[i] * numerator[i] / c_denominator[i] ;
            dn[i] = jnmx[i] * m * numerator[i] / d_denominator[i] ;
        }
    }

}

// -
