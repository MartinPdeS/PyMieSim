#pragma once

#include <experiment/polarization_set/polarization_set.h>
#include <experiment/base_set.h>
#include <experiment/material_set/material_set.h>

#include <single/detector/photodiode.h>
#include <single/detector/coherent_mode.h>
#include <single/detector/integrating_sphere.h>

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
};

class PhotodiodeSet : public BaseDetectorSet
{
    public:
        inline static const std::vector<std::string> attributes = {
            "sampling",
            "numerical_aperture",
            "cache_numerical_aperture",
            "phi_offset",
            "gamma_offset",
            "polarization_filter",
            "medium"
        };

        std::vector<unsigned> sampling;
        std::vector<double> numerical_aperture;
        std::vector<double> cache_numerical_aperture;
        std::vector<double> phi_offset;
        std::vector<double> gamma_offset;
        PolarizationSet polarization_filter_set;
        MediumSet medium;

        PhotodiodeSet() = default;

        PhotodiodeSet(
            const std::vector<unsigned> &sampling,
            const std::vector<double> &numerical_aperture,
            const std::vector<double> &cache_numerical_aperture,
            const std::vector<double> &phi_offset,
            const std::vector<double> &gamma_offset,
            const PolarizationSet &polarization_filter_set,
            const MediumSet &medium,
            const bool is_sequential
        );

        void update_shape() override;

        std::shared_ptr<BaseDetector> get_detector_by_index(long long flat_index) const override;

        std::shared_ptr<BaseDetector> get_detector_by_index_sequential(size_t index) const override;
};


class CoherentModeSet : public BaseDetectorSet
{
    public:
        inline static const std::vector<std::string> attributes = {
            "mode_numbers",
            "sampling",
            "numerical_aperture",
            "cache_numerical_aperture",
            "phi_offset",
            "gamma_offset",
            "polarization_filter",
            "rotation",
            "medium"
        };
        std::vector<std::string> mode_numbers;
        std::vector<unsigned> sampling;
        std::vector<double> numerical_aperture;
        std::vector<double> cache_numerical_aperture;
        std::vector<double> phi_offset;
        std::vector<double> gamma_offset;
        PolarizationSet polarization_filter_set;
        std::vector<double> rotation;
        MediumSet medium;
        bool coherent;
        bool mean_coupling;

        CoherentModeSet() = default;

        CoherentModeSet(
            const std::vector<std::string> &mode_numbers,
            const std::vector<unsigned> &sampling,
            const std::vector<double> &numerical_aperture,
            const std::vector<double> &cache_numerical_aperture,
            const std::vector<double> &phi_offset,
            const std::vector<double> &gamma_offset,
            const PolarizationSet &polarization_filter_set,
            const std::vector<double> &rotation,
            const MediumSet &medium,
            const bool &mean_coupling,
            const bool is_sequential
        );

        void update_shape() override;

        std::shared_ptr<BaseDetector> get_detector_by_index(long long flat_index) const override;

        std::shared_ptr<BaseDetector> get_detector_by_index_sequential(size_t index) const override;
};
