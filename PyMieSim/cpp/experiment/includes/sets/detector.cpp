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
            std::vector<double> cache_NA;
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
                const std::vector<double> &cache_NA,
                const std::vector<double> &phi_offset,
                const std::vector<double> &gamma_offset,
                const std::vector<double> &polarization_filter,
                const std::vector<double> &rotation,
                const bool &coherent,
                const bool &mean_coupling)
            : mode_numbers(mode_numbers), sampling(sampling), NA(NA), cache_NA(cache_NA), phi_offset(phi_offset), gamma_offset(gamma_offset),
              polarization_filter(polarization_filter), rotation(rotation), coherent(coherent), mean_coupling(mean_coupling)
              {
                update_shape();
                total_combinations = get_vector_sigma(shape);
              }

            void update_shape() override {
                this->shape = {
                    mode_numbers.size(),
                    sampling.size(),
                    NA.size(),
                    cache_NA.size(),
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
                    this->cache_NA[indices[3]],
                    this->phi_offset[indices[4]],
                    this->gamma_offset[indices[5]],
                    this->polarization_filter[indices[6]],
                    this->rotation[indices[7]],
                    this->coherent,
                    this->mean_coupling
                );

                detector.indices = indices;

                return detector;
            }

        void validate_sequential_data(const size_t expected_size) const {
            // Check each vector's size and throw an error with the specific vector name if sizes don't match
            if (this->mode_numbers.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. mode_numbers has a different size than expected size.");

            if (this->sampling.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. sampling has a different size than expected size.");

            if (this->NA.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. NA has a different size than expected size.");

            if (this->cache_NA.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. cache_NA has a different size than expected size.");

            if (this->phi_offset.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. phi_offset has a different size than expected size.");

            if (this->gamma_offset.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. gamma_offset has a different size than expected size.");

            if (this->polarization_filter.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. polarization_filter has a different size than expected size.");

            if (this->rotation.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. rotation has a different size than expected size.");
        }

            Detector get_detector_by_index_sequential(size_t index) const {

                return Detector(
                    this->mode_numbers[index],
                    this->sampling[index],
                    this->NA[index],
                    this->cache_NA[index],
                    this->phi_offset[index],
                    this->gamma_offset[index],
                    this->polarization_filter[index],
                    this->rotation[index],
                    this->coherent,
                    this->mean_coupling
                );
            }
    };
}
