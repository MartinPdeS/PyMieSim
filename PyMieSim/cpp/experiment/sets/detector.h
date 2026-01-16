#pragma once

#include "./base_set.h"
#include "single/detector/detector.h"


class DetectorSet : public BaseSet
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

        DetectorSet() = default;

        DetectorSet(const std::vector<std::string> &mode_numbers,
            const std::vector<unsigned> &sampling,
            const std::vector<double> &NA,
            const std::vector<double> &cache_NA,
            const std::vector<double> &phi_offset,
            const std::vector<double> &gamma_offset,
            const std::vector<double> &polarization_filter,
            const std::vector<double> &rotation,
            const bool &coherent,
            const bool &mean_coupling,
            const bool is_sequential
        )
        : BaseSet(is_sequential), mode_numbers(mode_numbers), sampling(sampling), NA(NA), cache_NA(cache_NA), phi_offset(phi_offset), gamma_offset(gamma_offset),
            polarization_filter(polarization_filter), rotation(rotation), coherent(coherent), mean_coupling(mean_coupling)
            {
            this->is_empty = false;
            this->update_shape();
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
            total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
        }

        Detector get_detector_by_index(long long flat_index) const {
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
            this->check_size(this->mode_numbers, expected_size, "mode_numbers");
            this->check_size(this->sampling, expected_size, "sampling");
            this->check_size(this->NA, expected_size, "NA");
            this->check_size(this->cache_NA, expected_size, "cache_NA");

            this->check_size(this->phi_offset, expected_size, "phi_offset");
            this->check_size(this->gamma_offset, expected_size, "gamma_offset");
            this->check_size(this->polarization_filter, expected_size, "polarization_filter");
            this->check_size(this->rotation, expected_size, "rotation");
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
