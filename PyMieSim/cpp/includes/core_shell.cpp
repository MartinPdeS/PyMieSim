#pragma once

#include "core_shell.h"

namespace CORESHELL
{

    void Scatterer::apply_medium(){
        this->core_index /= this->n_medium;
        this->shell_index /= this->n_medium;
        this->core_diameter *= this->n_medium;
        this->shell_width *= this->n_medium;
        this->shell_diameter *= this->n_medium;
    }

    void Scatterer::initialize(size_t max_order){
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

    void Scatterer::compute_max_order(size_t max_order){
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
        an.resize(max_order);
        bn.resize(max_order);

        complex128
            relative_index = this->shell_index / this->core_index,
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

        std::vector<complex128> pv(max_order + 1), pw(max_order + 1), py(max_order + 1), chv(max_order + 1), chw(max_order + 1), chy(max_order + 1), gsy(max_order + 1), gs1y(max_order + 1);
        std::vector<complex128> p1y(max_order + 2), ch1y(max_order + 2);

        p1y[0] = sin(this->x_shell);
        ch1y[0] = cos(this->x_shell);

        for (size_t order = 0; order < max_order + 1; order++){
            double nu = order + 1.5 ;
            pw[order] = sw * compute_Jn(nu, w);
            pv[order] = sv * compute_Jn(nu, v);
            py[order] = sy * compute_Jn(nu, this->x_shell);

            chv[order] = -sv * compute_Yn(nu, v);
            chw[order] = -sw * compute_Yn(nu, w);
            chy[order] = -sy * compute_Yn(nu, this->x_shell);

            p1y[order + 1] = py[order];
            ch1y[order + 1] = chy[order];
            gsy[order] = py[order]  - JJ * chy[order];
            gs1y[order] = p1y[order] - JJ * ch1y[order];
        }

        std::vector<complex128> Du(nmx, 0.0), Dv(nmx, 0.0), Dw(nmx, 0.0);

        for (int i = nmx-1; i > 1; i--){
            Du[i-1] = (double)i / u - 1.0 / (Du[i] + (double)i / u);
            Dv[i-1] = (double)i / v - 1.0 / (Dv[i] + (double)i / v);
            Dw[i-1] = (double)i / w - 1.0 / (Dw[i] + (double)i / w);
        }

        Du.erase(Du.begin());
        Dv.erase(Dv.begin());
        Dw.erase(Dw.begin());

        std::vector<complex128> uu(max_order), vv(max_order), fv(max_order), dns(max_order), gns(max_order), a1(max_order), b1(max_order);

        for (size_t order=0; order < max_order; order++){
            double idx = order + 1;
            uu[order] = relative_index * Du[order] - Dv[order];
            vv[order] = Du[order] / relative_index - Dv[order];
            fv[order] = pv[order] / chv[order];
            dns[order] = ((uu[order] * fv[order] / pw[order]) / (uu[order] * (pw[order] - chw[order] * fv[order]) + pw[order] / pv[order] / chv[order])) + Dw[order];
            gns[order] = ((vv[order] * fv[order] / pw[order]) / (vv[order] * (pw[order] - chw[order] * fv[order]) + pw[order] / pv[order] / chv[order])) + Dw[order];
            a1[order] = dns[order] / shell_index + idx / x_shell;
            b1[order] = shell_index * gns[order] + idx / x_shell;
            an[order] = (py[order] * a1[order] - p1y[order]) / (gsy[order] * a1[order] - gs1y[order]);
            bn[order] = (py[order] * b1[order] - p1y[order]) / (gsy[order] * b1[order] - gs1y[order]);
        }
    }

}

// -
