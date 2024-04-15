#pragma once

#include "core_shell.h"

namespace CORESHELL
{

    void Scatterer::initialize(size_t &max_order){
        this->apply_medium();
        compute_size_parameter();
        compute_max_order(max_order);
        compute_area();
        compute_an_bn();
    }

    void Scatterer::compute_size_parameter(){
        size_parameter = source.k * this->shell_diameter / 2;
        this->x_shell = source.k * this->shell_diameter / 2.0;
        this->x_core = source.k * this->core_diameter / 2.0;
    }

    void Scatterer::compute_max_order(size_t &max_order){
        if (max_order == 0)
            this->max_order = (size_t) (2 + size_parameter + 4 * pow(size_parameter, 1./3.));
        else
            this->max_order = max_order;
    }

    void Scatterer::compute_area(){

        this->area = PI * pow(this->shell_diameter/2.0, 2);
    }

    void Scatterer::compute_an_bn()
    {
        an = std::vector<complex128>(max_order);
        bn = std::vector<complex128>(max_order);

        complex128
            m = this->shell_index / this->core_index,
            u = this->core_index * this->x_core,
            v = this->shell_index * this->x_core,
            w = this->shell_index * this->x_shell;

        complex128
            sv = sqrt(0.5 * PI * v),
            sw = sqrt(0.5 * PI * w),
            sy = sqrt(0.5 * PI * this->x_shell);

        size_t
            mx = (size_t) std::max( abs( this->core_index * this->x_shell ), abs( this->shell_index*this->x_shell ) ),
            nmx  = (size_t) ( std::max( max_order, mx ) + 16. )  ;

        std::vector<complex128> pv, pw, py, chv, chw, chy, p1y, ch1y, gsy, gs1y;

        p1y.push_back( sin( this->x_shell ) ) ;
        ch1y.push_back( cos( this->x_shell ) ) ;

        for (size_t i=0; i<max_order+1; i++){
            double nu = i + 1.5 ;
            pw.push_back( sw * compute_Jn(nu, w) );
            pv.push_back( sv * compute_Jn(nu, v) );
            py.push_back( sy * compute_Jn(nu, this->x_shell) );

            chv.push_back( -sv * compute_Yn(nu, v) );
            chw.push_back( -sw * compute_Yn(nu, w) );
            chy.push_back( -sy * compute_Yn(nu, this->x_shell) );

            p1y.push_back ( py[i]  );
            ch1y.push_back( chy[i] );
            gsy.push_back ( py[i]  - JJ * chy[i]  );
            gs1y.push_back ( p1y[i] - JJ * ch1y[i] );
        }

        std::vector<complex128>
            Du = std::vector<complex128>(nmx, 0.),
            Dv = std::vector<complex128>(nmx, 0.),
            Dw = std::vector<complex128>(nmx, 0.);

        for (int i = nmx-1; i > 1; i--){
            Du[i-1] = (double)i / u -1. / (Du[i] + (double)i / u);
            Dv[i-1] = (double)i / v -1. / (Dv[i] + (double)i / v);
            Dw[i-1] = (double)i / w -1. / (Dw[i] + (double)i / w);
        }

        Du.erase(Du.begin());
        Dv.erase(Dv.begin());
        Dw.erase(Dw.begin());

        std::vector<complex128> uu, vv, fv, dns, gns, a1, b1;
        for (size_t i=0; i<max_order; i++){
            double n = (double) (i+1);
            uu.push_back ( m * Du[i] - Dv[i]  );
            vv.push_back ( Du[i] / m - Dv[i] );
            fv.push_back ( pv[i] / chv[i]    );
            dns.push_back( ( ( uu[i] * fv[i] / pw[i] ) / ( uu[i] * ( pw[i] - chw[i] * fv[i] ) + ( pw[i] / pv[i] ) / chv[i] ) ) + Dw[i] );
            gns.push_back( ( ( vv[i] * fv[i] / pw[i] ) / ( vv[i] * ( pw[i] - chw[i] * fv[i] ) + ( pw[i] / pv[i] ) / chv[i] ) ) + Dw[i] );
            a1.push_back ( dns[i] / this->shell_index + n / this->x_shell );
            b1.push_back ( this->shell_index * gns[i] + n / this->x_shell );
            an[i] = ( py[i] * a1[i] - p1y[i] ) / ( gsy[i] * a1[i] - gs1y[i] ) ;
            bn[i] = ( py[i] * b1[i] - p1y[i] ) / ( gsy[i] * b1[i] - gs1y[i] ) ;
        }
    }


    std::tuple<std::vector<complex128>, std::vector<complex128>> Scatterer::compute_s1s2(const std::vector<double> &phi) const {
        std::vector<complex128>
            S1(phi.size(), 0.0),
            S2(phi.size(), 0.0);

        std::vector<double> prefactor = get_prefactor();

        std::vector<double> mu;
        mu.reserve(phi.size());

        for (double phi : phi)
            mu.push_back( cos( phi-PI / 2.0 ) );


        for (unsigned int i = 0; i < phi.size(); i++){
            auto [pin, taun] = VSH::SPHERICAL::MiePiTau(mu[i], max_order);

            for (unsigned int m = 0; m < max_order ; m++){
                S1[i] += prefactor[m] * ( this->an[m] * pin[m] +  this->bn[m] * taun[m] );
                S2[i] += prefactor[m] * ( this->an[m] * taun[m] + this->bn[m] * pin[m]  );
            }
        }

        return std::make_tuple(S1, S2);
    }


    double Scatterer::get_Qsca() const {
        double value = 0;

        for(size_t it = 0; it < max_order; ++it){
            double n = (double) it + 1;
            value += (2.* n + 1.) * ( pow( std::abs(this->an[it]), 2) + pow( std::abs(this->bn[it]), 2)  );
        }

        return value * 2. / pow( size_parameter, 2.);
    }


    double Scatterer::get_Qext() const {
        double value = 0;

        for(size_t it = 0; it < max_order; ++it){
            double n = (double) it + 1;
            value += (2.* n + 1.) * std::real( this->an[it] + this->bn[it] );
        }
        return value * 2. / pow( size_parameter, 2.);
    }


    double Scatterer::get_Qback() const {
        complex128 value = 0;

        for(size_t it = 0; it < max_order-1; ++it){
            double n = (double) it + 1;
            value += (2. * n + 1) * pow(-1., n) * ( this->an[it] - this->bn[it] ) ;
        }

        value *= pow( std::abs(value), 2. ) / pow( size_parameter, 2. );
        return std::abs(value);
    }

    double Scatterer::get_g() const {
        double Qsca = get_Qsca(),
        value = 0;

        for(size_t it = 0; it < max_order-1; ++it){
            double n = (double) it + 1;
            value += ( n * (n + 2.) / (n + 1.) ) * std::real(this->an[it] * std::conj(this->an[it+1]) + this->bn[it] * std::conj(this->bn[it+1]) );
            value += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( this->an[it] * std::conj(this->bn[it]) );
        }

        return value * 4. / ( Qsca * pow(size_parameter, 2) );
    }

}

// -
