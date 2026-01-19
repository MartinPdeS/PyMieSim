#pragma once

#include "./base_set.h"
#include "single/detector/detector.h"

class BaseDetectorSet : public BaseSet
{
    public:
        virtual ~BaseDetectorSet() = default;
        /**
         * @brief Retrieves a detector configuration by its flat index.
         * @param flat_index The flat index of the desired detector configuration.
         * @return The corresponding BaseDetector configuration.
         */
        virtual std::shared_ptr<BaseDetector> get_detector_by_index(long long flat_index) const = 0;

        /**
         * @brief Retrieves a detector configuration by its index in sequential mode.
         * @param index The index of the desired detector configuration.
         * @return The corresponding Photodiode detector configuration.
         */
        virtual std::shared_ptr<BaseDetector> get_detector_by_index_sequential(size_t index) const = 0;

        /**
         * @brief Validates the sizes of the sequential data vectors.
         * @param expected_size The expected size for each vector.
         */
        virtual void validate_sequential_data(const size_t expected_size) const = 0;
};

class PhotodiodeSet : public BaseDetectorSet
{
    public:
        std::vector<unsigned> sampling;
        std::vector<double> NA;
        std::vector<double> cache_NA;
        std::vector<double> phi_offset;
        std::vector<double> gamma_offset;
        std::vector<double> polarization_filter;
        std::vector<double> medium_refractive_index;

        PhotodiodeSet() = default;

        PhotodiodeSet(
            const std::vector<unsigned> &sampling,
            const std::vector<double> &NA,
            const std::vector<double> &cache_NA,
            const std::vector<double> &phi_offset,
            const std::vector<double> &gamma_offset,
            const std::vector<double> &polarization_filter,
            const std::vector<double> &medium_refractive_index,
            const bool is_sequential
        )
        :   sampling(sampling),
            NA(NA),
            cache_NA(cache_NA),
            phi_offset(phi_offset),
            gamma_offset(gamma_offset),
            polarization_filter(polarization_filter),
            medium_refractive_index(medium_refractive_index)
            {
                this->is_sequential = is_sequential;
                this->is_empty = false;
                this->update_shape();
            }

        void update_shape() override {
            this->shape = {
                sampling.size(),
                NA.size(),
                cache_NA.size(),
                phi_offset.size(),
                gamma_offset.size(),
                polarization_filter.size(),
                medium_refractive_index.size()
            };
            total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
        }

        std::shared_ptr<BaseDetector> get_detector_by_index(long long flat_index) const override{
            std::vector<size_t> indices = this->calculate_indices(flat_index);

            std::shared_ptr<Photodiode> detector = std::make_shared<Photodiode>(
                this->sampling[indices[0]],
                this->NA[indices[1]],
                this->cache_NA[indices[2]],
                this->phi_offset[indices[3]],
                this->gamma_offset[indices[4]],
                this->polarization_filter[indices[5]],
                this->medium_refractive_index[indices[6]]
            );

            detector->indices = indices;

            return detector;
        }

        void validate_sequential_data(const size_t expected_size) const override {
            // Check each vector's size and throw an error with the specific vector name if sizes don't match
            this->check_size(this->sampling, expected_size, "sampling");
            this->check_size(this->NA, expected_size, "NA");
            this->check_size(this->cache_NA, expected_size, "cache_NA");

            this->check_size(this->phi_offset, expected_size, "phi_offset");
            this->check_size(this->gamma_offset, expected_size, "gamma_offset");
            this->check_size(this->polarization_filter, expected_size, "polarization_filter");
            this->check_size(this->medium_refractive_index, expected_size, "medium_refractive_index");
        }

        std::shared_ptr<BaseDetector> get_detector_by_index_sequential(size_t index) const override {
            return std::make_shared<Photodiode>(
                this->sampling[index],
                this->NA[index],
                this->cache_NA[index],
                this->phi_offset[index],
                this->gamma_offset[index],
                this->polarization_filter[index],
                this->medium_refractive_index[index]
            );
        }
};


class CoherentModeSet : public BaseDetectorSet
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
        std::vector<double> medium_refractive_index;
        bool coherent;
        bool mean_coupling;

        CoherentModeSet() = default;

        CoherentModeSet(
            const std::vector<std::string> &mode_numbers,
            const std::vector<unsigned> &sampling,
            const std::vector<double> &NA,
            const std::vector<double> &cache_NA,
            const std::vector<double> &phi_offset,
            const std::vector<double> &gamma_offset,
            const std::vector<double> &polarization_filter,
            const std::vector<double> &rotation,
            const std::vector<double> &medium_refractive_index,
            const bool &mean_coupling,
            const bool is_sequential
        )
        :   mode_numbers(mode_numbers),
            sampling(sampling),
            NA(NA),
            cache_NA(cache_NA),
            phi_offset(phi_offset),
            gamma_offset(gamma_offset),
            polarization_filter(polarization_filter),
            rotation(rotation),
            medium_refractive_index(medium_refractive_index),
            mean_coupling(mean_coupling)
            {
                this->is_sequential = is_sequential;
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
                rotation.size(),
                medium_refractive_index.size()
            };
            total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
        }

        std::shared_ptr<BaseDetector> get_detector_by_index(long long flat_index) const override {
            std::vector<size_t> indices = this->calculate_indices(flat_index);
            std::shared_ptr<CoherentMode> detector = std::make_shared<CoherentMode>(
                this->mode_numbers[indices[0]],
                this->sampling[indices[1]],
                this->NA[indices[2]],
                this->cache_NA[indices[3]],
                this->phi_offset[indices[4]],
                this->gamma_offset[indices[5]],
                this->polarization_filter[indices[6]],
                this->rotation[indices[7]],
                this->mean_coupling,
                this->medium_refractive_index[indices[8]]
            );

            detector->indices = indices;

            return detector;
        }

        void validate_sequential_data(const size_t expected_size) const override {
            // Check each vector's size and throw an error with the specific vector name if sizes don't match
            this->check_size(this->mode_numbers, expected_size, "mode_numbers");
            this->check_size(this->sampling, expected_size, "sampling");
            this->check_size(this->NA, expected_size, "NA");
            this->check_size(this->cache_NA, expected_size, "cache_NA");

            this->check_size(this->phi_offset, expected_size, "phi_offset");
            this->check_size(this->gamma_offset, expected_size, "gamma_offset");
            this->check_size(this->polarization_filter, expected_size, "polarization_filter");
            this->check_size(this->rotation, expected_size, "rotation");
            this->check_size(this->medium_refractive_index, expected_size, "medium_refractive_index");
        }

        std::shared_ptr<BaseDetector> get_detector_by_index_sequential(size_t index) const override {
            return std::make_shared<CoherentMode>(
                this->mode_numbers[index],
                this->sampling[index],
                this->NA[index],
                this->cache_NA[index],
                this->phi_offset[index],
                this->gamma_offset[index],
                this->polarization_filter[index],
                this->rotation[index],
                this->mean_coupling,
                this->medium_refractive_index[index]
            );
        }
};
