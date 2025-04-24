#pragma once

#include <complex>
#include "scatterer/base_scatterer/base_scatterer.h"


using complex128 = std::complex<double>;


class Sphere: public BaseScatterer
{
    public:
        double diameter;
        complex128 refractive_index;

        Sphere(const double diameter, const complex128 refractive_index, const double medium_refractive_index, const BaseSource &source, size_t max_order = 0, bool compute_cn_dn = false)
        : BaseScatterer(max_order, source, medium_refractive_index), diameter(diameter), refractive_index(refractive_index)
        {
            this->compute_area();
            this->compute_size_parameter();
            this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
            this->compute_an_bn(this->max_order);

            if (compute_cn_dn)
                this->compute_cn_dn(this->max_order);
        }

        void compute_size_parameter() override {
            this->size_parameter = source.wavenumber * this->diameter / 2 * this->medium_refractive_index;
            this->size_parameter_squared = pow(this->size_parameter, 2);
        }

        void compute_area() override {
            this->area = PI * std::pow(this->diameter / 2.0, 2);
        }

        void compute_an_bn(const size_t max_order);
        void compute_cn_dn(const size_t max_order);

        double get_Qsca() const override;
        double get_Qext() const override;
        double get_Qback() const override;
        double get_g() const override;

        std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const override {
            std::vector<complex128> S1, S2;

            S1.reserve(phi.size());
            S2.reserve(phi.size());

            std::vector<double> mu, prefactor = get_prefactor();

            mu.reserve(phi.size());

            for (const double phi : phi)
                mu.push_back( cos( phi - PI / 2.0 ) );

            for (size_t i = 0; i < phi.size(); i++){
                auto [pin, taun] = this->get_pi_tau(mu[i], max_order);
                complex128 S1_temp = 0., S2_temp = 0.;

                for (size_t m = 0; m < max_order ; m++){
                    S1_temp += prefactor[m] * ( this->an[m] * pin[m] +  this->bn[m] * taun[m] );
                    S2_temp += prefactor[m] * ( this->an[m] * taun[m] + this->bn[m] * pin[m]  );
                }
                S1.push_back(S1_temp);
                S2.push_back(S2_temp);
            }

            return std::make_tuple(std::move(S1), std::move(S2));
        }

};

// -
