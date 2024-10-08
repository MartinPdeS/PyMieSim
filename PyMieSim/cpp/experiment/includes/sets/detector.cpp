#pragma once


#include "experiment/includes/sets/base.cpp"
#include "single/includes/detectors.cpp"


using complex128 = std::complex<double>;

namespace DETECTOR
{
    class Set : public BaseSet
    {
        public:
            std::vector<std::string> mode_numbers;
            std::vector<unsigned> sampling;
            std::vector<double> NA;
            std::vector<double> phi_offset;
            std::vector<double> gamma_offset;
            std::vector<double> polarization_filter;
            std::vector<double> rotation;
            bool coherent;
            bool mean_coupling;

            Set() = default;

            Set(const std::vector<std::string> &mode_numbers,
                const std::vector<unsigned> &sampling,
                const std::vector<double> &NA,
                const std::vector<double> &phi_offset,
                const std::vector<double> &gamma_offset,
                const std::vector<double> &polarization_filter,
                const std::vector<double> &rotation,
                const bool &coherent,
                const bool &mean_coupling)
            : mode_numbers(mode_numbers), sampling(sampling), NA(NA), phi_offset(phi_offset), gamma_offset(gamma_offset),
              polarization_filter(polarization_filter), rotation(rotation), coherent(coherent), mean_coupling(mean_coupling)
              {
                update_shape();
                total_combinations = get_vector_sigma(shape);
              }

            void update_shape() override {
                shape = {
                    mode_numbers.size(),
                    sampling.size(),
                    NA.size(),
                    phi_offset.size(),
                    gamma_offset.size(),
                    polarization_filter.size(),
                    rotation.size()
                };
            }

            Detector get_detector_by_index(size_t flat_index) const {

                std::vector<size_t> indices = this->calculate_indices(flat_index);

                Detector detector(
                    this->mode_numbers[indices[0]],
                    this->sampling[indices[1]],
                    this->NA[indices[2]],
                    this->phi_offset[indices[3]],
                    this->gamma_offset[indices[4]],
                    this->polarization_filter[indices[5]],
                    this->rotation[indices[6]],
                    this->coherent,
                    this->mean_coupling
                );

                detector.indices = indices;

                return detector;
            }
    };
}
