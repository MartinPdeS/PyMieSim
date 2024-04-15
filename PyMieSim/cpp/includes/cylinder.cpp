#pragma once

#include "cylinder.h"

namespace CYLINDER
{
    void Scatterer::initialize(size_t &max_order){
        // this->apply_medium();
        this->compute_size_parameter();
        this->compute_max_order(max_order);
        this->compute_area();
        this->compute_an_bn();

    }

    void Scatterer::compute_max_order(size_t &max_order){
        if (max_order == 0)
            this->max_order  = (size_t) (2 + this->size_parameter + 4 * pow(this->size_parameter, 1./3.)) + 16;
        else
            this->max_order = max_order;
    }

    void Scatterer::compute_size_parameter(){
        this->size_parameter = PI * this->diameter / source.wavelength;
    }

    void Scatterer::compute_area(){
        this->area = this->diameter;
    }

    double Scatterer::get_g() const {
        return get_g_with_fields(1000, 1.0);
    }

    double Scatterer::get_Qsca() const {
        complex128 Qsca1=0, Qsca2=0;

        for(size_t it = 1; it < max_order; it++)
        {
            Qsca1 +=  pow( std::abs(this->a1n[it]), 2 ) + pow( std::abs(this->b1n[it]), 2 ) ;
            Qsca2 +=  pow( std::abs(this->a2n[it]), 2 ) + pow( std::abs(this->b2n[it]), 2 ) ;
        }

        Qsca1 =  2. / size_parameter * ( 2.0 * Qsca1 + pow( abs(this->b1n[0]), 2 ) );
        Qsca2 =  2. / size_parameter * ( 2.0 * Qsca2 + pow( abs(this->a2n[0]), 2 ) );

        return process_polarization(Qsca1, Qsca2);
    }

    double Scatterer::get_Qext() const {
        complex128 Qext1 = 0, Qext2 = 0;

        for(size_t it = 1; it < max_order; ++it){
            Qext1 += this->b1n[it];
            Qext2 += this->a2n[it];
        }

        Qext1 = 2. / size_parameter * std::real( this->b1n[0] + 2.0 * Qext1 );
        Qext2 = 2. / size_parameter * std::real( this->a1n[0] + 2.0 * Qext2 );

        return this->process_polarization(Qext1, Qext2);
    }

    double Scatterer::process_polarization(complex128 &value_0, complex128& value_1) const {
        // if (source.is_polarized == false)
        //     return 0.5 * ( abs( value1 ) + abs( Value0 ) );
        // else
            return abs( value_1 ) * pow(abs(source.jones_vector[0]), 2) + abs( value_0 ) * pow(abs(source.jones_vector[1]), 2);
    }

    void Scatterer::compute_an_bn() {
        this->a1n = std::vector<complex128>(max_order);
        this->b1n = std::vector<complex128>(max_order);

        this->a2n = std::vector<complex128>(max_order);
        this->b2n = std::vector<complex128>(max_order);

        double x = size_parameter;

        complex128
            m = this->index / this->n_medium,
            z = m * x,
            numerator,
            denominator;

        std::vector<complex128>
            J_z(max_order+1),
            J_z_p(max_order+1),
            J_x(max_order+1),
            J_x_p(max_order+1),
            H_x(max_order+1),
            H_x_p(max_order+1);

        for (size_t n = 0; n < max_order+1; ++n){
            double nd = (double) n ;

            J_z[n] = compute_Jn(nd, z);
            J_z_p[n] = compute_Jn_p(nd, z);
            J_x[n] = compute_Jn(nd, x);
            J_x_p[n] = compute_Jn_p(nd, x);
            H_x[n] = compute_H1(nd, x);
            H_x_p[n] = compute_H1_p(nd, x);
        }

        for (size_t n = 0; n < (size_t) max_order; n++){
            numerator = m * J_z[n] * J_x_p[n] - J_z_p[n] * J_x[n];
            denominator = m * J_z[n] * H_x_p[n] - J_z_p[n] * H_x[n];
            this->a1n[n] = 0.0 ;
            this->a2n[n] = numerator/denominator ;

            numerator = J_z[n] * J_x_p[n] - m * J_z_p[n] * J_x[n];
            denominator = J_z[n] * H_x_p[n] - m * J_z_p[n] * H_x[n];
            this->b1n[n] = numerator/denominator ;
            this->b2n[n] = 0.0 ;
        }
    }

    std::tuple<std::vector<complex128>, std::vector<complex128>> Scatterer::compute_s1s2(const std::vector<double> &phi) const{
        std::vector<complex128>
        T1(phi.size()),
        T2(phi.size());

        for (unsigned int i = 0; i < phi.size(); i++){
            T1[i] = this->b1n[0];
            T2[i] = this->a2n[0];
            for (size_t order = 1; order < max_order ; order++){
                T1[i] += 2.0 * this->b1n[order] * cos(order * (PI - (phi[i] + PI/2.0) ) );
                T2[i] += 2.0 * this->a2n[order] * cos(order * (PI - (phi[i] + PI/2.0) ) );
            }
        }

        return std::make_tuple( T1, T2 )  ;
    }

}
