#pragma once

#include "experiment/sets/base_set.h"
#include "single/source/source.h"
#include "./polarization_set.h"

using complex128 = std::complex<double>;


class BaseSourceSet : public BaseSet {  // Abstract Class for Sources
public:
    std::vector<double> wavelength;

    std::shared_ptr<PolarizationSet> polarization;
    std::vector<double> numerical_aperture;
    std::vector<double> optical_power;
    std::vector<double> amplitude;

    BaseSourceSet() = default;

    virtual ~BaseSourceSet() = default;

    virtual std::shared_ptr<BaseSource> get_source_by_index(const size_t) const = 0;

    virtual std::shared_ptr<BaseSource> get_source_by_index_sequential(const size_t) const = 0;

    virtual void validate_sequential_data(const size_t) const = 0;
};


class GaussianSourceSet : public BaseSourceSet
{
public:
    std::vector<std::string> attributes = {
        "wavelength",
        "polarization",
        "numerical_aperture",
        "optical_power"
    };
    GaussianSourceSet() = default;

    GaussianSourceSet(
        const std::vector<double>& wavelength,
        const std::shared_ptr<PolarizationSet>& polarization,
        const std::vector<double>& numerical_aperture,
        const std::vector<double>& optical_power,
        const bool is_sequential)
    {
        this->wavelength = wavelength;
        this->polarization = polarization;
        this->numerical_aperture = numerical_aperture;
        this->optical_power = optical_power;
        this->is_sequential = is_sequential;

        this->shape = {wavelength.size(), polarization->number_of_states(), numerical_aperture.size(), optical_power.size()};
        this->total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
    }


    std::shared_ptr<BaseSource> get_source_by_index(const size_t flat_index) const override;

    std::shared_ptr<BaseSource> get_source_by_index_sequential(const size_t index) const override;

    void validate_sequential_data(const size_t expected_size) const override;

};


class PlaneWaveSourceSet : public BaseSourceSet
{
public:
    std::vector<std::string> attributes = {
        "wavelength",
        "polarization",
        "amplitude"
    };
    PlaneWaveSourceSet() = default;

    PlaneWaveSourceSet(const std::vector<double>& wavelength,
        const std::shared_ptr<PolarizationSet>& polarization,
        const std::vector<double>& amplitude,
        const bool is_sequential)
    {
        this->wavelength = wavelength;
        this->polarization = polarization;
        this->amplitude = amplitude;
        this->is_sequential = is_sequential;

        this->shape = {wavelength.size(), polarization->number_of_states(), amplitude.size()};
        this->total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
    }

    std::shared_ptr<BaseSource> get_source_by_index(const size_t flat_index) const override;

    std::shared_ptr<BaseSource> get_source_by_index_sequential(const size_t index) const override;

    void validate_sequential_data(const size_t expected_size) const override;
};
