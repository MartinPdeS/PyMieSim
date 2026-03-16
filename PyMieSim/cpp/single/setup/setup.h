#pragma once

#define RETURN_IF_MATCH(name) \
    if (data_name == #name) return scatterer->get_##name();

#include <functional>
#include <memory>
#include <cstdarg>
#include <cstdio>
#include <single/scatterer/sphere/sphere.h>
#include <single/scatterer/cylinder/cylinder.h>
#include <single/scatterer/coreshell/coreshell.h>
#include <single/source/source.h>

#include <single/detector/photodiode.h>
#include <single/detector/coherent_mode.h>
#include <single/detector/integrating_sphere.h>

#include <utils/constants.h>

class Setup
{
    public:
        const std::shared_ptr<BaseSource> source;
        std::shared_ptr<BaseScatterer> scatterer;
        const std::shared_ptr<BaseDetector> detector;

        bool debug_mode = false;

        Setup() = default;
        Setup(
            std::shared_ptr<BaseScatterer> _scatterer,
            std::shared_ptr<BaseSource> _source,
            std::shared_ptr<BaseDetector> _detector = nullptr,
            bool _debug_mode = false)
            : source(_source), scatterer(_scatterer), detector(_detector), debug_mode(_debug_mode)
        {
            scatterer->init(source, 0);  // max_order = 0 for automatic selection
            if (detector) {
                detector->initialize_mesh(scatterer);
                detector->medium->initialize(source->wavelength);
            }

        }

        double get(const std::string& data_name) const {
            if (data_name == "coupling") {
                if (!detector) {
                    throw std::logic_error("Detector is not defined in the setup. Coupling cannot be computed.");
                }
                return detector->get_coupling(scatterer, source);

            }

            RETURN_IF_MATCH(cross_section)

            RETURN_IF_MATCH(size_parameter)

            RETURN_IF_MATCH(Qsca)
            RETURN_IF_MATCH(Csca)

            RETURN_IF_MATCH(Qext)
            RETURN_IF_MATCH(Cext)

            RETURN_IF_MATCH(Qabs)
            RETURN_IF_MATCH(Cabs)

            RETURN_IF_MATCH(Qback)
            RETURN_IF_MATCH(Cback)

            RETURN_IF_MATCH(Qforward)
            RETURN_IF_MATCH(Cforward)

            RETURN_IF_MATCH(Qratio)
            RETURN_IF_MATCH(Cratio)

            RETURN_IF_MATCH(Qpr)
            RETURN_IF_MATCH(Cpr)

            RETURN_IF_MATCH(g)
            if (data_name == "g_with_farfields") return scatterer->get_g_with_farfields(this->source, 1000);

            throw std::invalid_argument("Unknown data name: " + data_name);
        }

        std::vector<double> get_poynting_field(double distance) const
        {
            if (!detector) {
                throw std::logic_error("Detector is not defined in the setup. Poynting field cannot be computed.");
            }

            auto [theta_field, phi_field] = this->scatterer->get_unstructured_farfields(
                this->detector->fibonacci_mesh,
                distance,
                this->source
            );

            std::vector<double> poynting(phi_field.size());

            for (size_t i = 0; i < phi_field.size(); ++i) {
                const double fields_squared = std::norm(phi_field[i]) + std::norm(theta_field[i]);
                poynting[i] = (
                    this->scatterer->medium->get_refractive_index() *
                    Constants::EPSILON0 *
                    Constants::LIGHT_SPEED *
                    fields_squared
                );
            }

            return poynting;
        }

    //------------------------ S1S2 --------------------------
    std::pair<std::vector<complex128>, std::vector<complex128>>
    get_s1s2(const std::vector<double>& angles) const
    {
        std::vector<double> phi_values = angles;
        for (size_t i = 0; i < phi_values.size(); ++i)
            phi_values[i] += Constants::PI / 2.0;

        return scatterer->compute_s1s2(phi_values);
    }

    //------------------------ FARFIELDS --------------------------
    std::pair<std::vector<complex128>, std::vector<complex128>>
    get_unstructured_farfields(
        const std::vector<double>& phi,
        const std::vector<double>& theta,
        double distance
    ) const {
        return this->scatterer->get_unstructured_farfields(
            phi,
            theta,
            distance,
            this->source
        );
    }

    std::tuple<
        std::vector<complex128>,
        std::vector<complex128>,
        FullSteradian
    >
    get_structured_farfields(
        const size_t sampling,
        double distance
    ) const {
        return this->scatterer->get_structured_farfields(
            sampling,
            distance,
            this->source
        );
    }

    //------------------------ NEARFIELDS --------------------------
    std::vector<complex128>
    get_scattered_nearfields(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
        const std::string& field_type
    ) const {
        return this->scatterer->get_scattered_nearfields(
            x,
            y,
            z,
            field_type,
            this->source
        );
    }

    std::vector<complex128>
    get_incident_nearfields(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
        const std::string& field_type
    ) const {
        return this->scatterer->compute_incident_nearfields(
            x,
            y,
            z,
            field_type,
            this->source
        );
    }


    std::vector<complex128> get_total_nearfields(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
        const std::string& field_type
    ) const {
        return this->scatterer->get_total_nearfields(
            x,
            y,
            z,
            field_type,
            this->source
        );
    }

    //------------------------ STOKES --------------------------
    std::tuple<
        std::vector<double>,
        std::vector<double>,
        std::vector<double>,
        std::vector<double>
    >
    get_unstructured_stokes(
        const std::vector<double>& phi,
        const std::vector<double>& theta,
        const double distance
    ) const {
        return this->scatterer->get_unstructured_stokes(
            phi,
            theta,
            distance,
            this->source
        );
    }

    std::tuple<
        std::vector<double>,
        std::vector<double>,
        std::vector<double>,
        std::vector<double>,
        FullSteradian
    >
    get_structured_stokes(
        const size_t sampling,
        const double distance
    ) const {
        return this->scatterer->get_structured_stokes(
            sampling,
            distance,
            this->source
        );
    }



    //------------------------ SPF --------------------------
    std::pair<std::vector<double>, FullSteradian>
    get_structured_spf(
        const size_t sampling,
        const double distance = 1.0
    ) const {
        return this->scatterer->get_structured_spf(
            this->source,
            sampling,
            distance
        );
    }

    std::vector<double>
    get_unstructured_spf(
        const std::vector<double>& phi,
        const std::vector<double>& theta,
        const double distance = 1.0
    ) const {
        return this->scatterer->get_unstructured_spf(
            this->source,
            phi,
            theta,
            distance
        );
    }


};



#undef RETURN_IF_MATCH