#pragma once

#include "base_set.h"
#include "single/source/source.h"

using complex128 = std::complex<double>;


class BaseSourceSet : public BaseSet{  // Abstract Class for Sources
public:
    std::vector<double> wavelength;
    std::vector<std::vector<complex128>> jones_vector;
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
    GaussianSourceSet() = default;

    GaussianSourceSet(
        const std::vector<double>& wavelength,
        const std::vector<std::vector<complex128>>& jones_vector,
        const std::vector<double>& numerical_aperture,
        const std::vector<double>& optical_power,
        const bool is_sequential)
    {
        this->wavelength = wavelength;
        this->jones_vector = jones_vector;
        this->numerical_aperture = numerical_aperture;
        this->optical_power = optical_power;
        this->is_sequential = is_sequential;

        this->shape = {wavelength.size(), jones_vector.size(), numerical_aperture.size(), optical_power.size()};
        this->total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
    }


    std::shared_ptr<BaseSource> get_source_by_index(const size_t flat_index) const override {
        std::vector<size_t> indices = calculate_indices(flat_index);

        std::shared_ptr<BaseSource> source = std::make_shared<Gaussian>(
            this->wavelength[indices[0]],
            this->jones_vector[indices[1]],
            this->numerical_aperture[indices[2]],
            this->optical_power[indices[3]]
        );

        source->indices = indices;
        source->wavelength_index = indices[0];

        return source;
    }

    std::shared_ptr<BaseSource> get_source_by_index_sequential(const size_t index) const override {
        std::shared_ptr<BaseSource> source = std::make_shared<Gaussian>(
            this->wavelength[index],
            this->jones_vector[index],
            this->numerical_aperture[index],
            this->optical_power[index]
        );

        return source;
    }

    void validate_sequential_data(const size_t expected_size) const override {
        // Check each vector's size and throw an error with the specific vector name if sizes don't match
        this->check_size(this->wavelength, expected_size, "wavelength");
        this->check_size(this->jones_vector, expected_size, "jones_vector");
        this->check_size(this->numerical_aperture, expected_size, "numerical_aperture");
        this->check_size(this->optical_power, expected_size, "optical_power");
    }

};


class PlaneWaveSourceSet : public BaseSourceSet
{
public:
    PlaneWaveSourceSet() = default;

    PlaneWaveSourceSet(const std::vector<double>& wavelength,
        const std::vector<std::vector<complex128>>& jones_vector,
        const std::vector<double>& amplitude,
        const bool is_sequential)
    {
        this->wavelength = wavelength;
        this->jones_vector = jones_vector;
        this->amplitude = amplitude;
        this->is_sequential = is_sequential;

        this->shape = {wavelength.size(), jones_vector.size(), amplitude.size()};
        this->total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
    }

    std::shared_ptr<BaseSource> get_source_by_index(const size_t flat_index) const override {
        std::vector<size_t> indices = calculate_indices(flat_index);

        std::shared_ptr<BaseSource> source = std::make_shared<Planewave>(
            this->wavelength[indices[0]],
            this->jones_vector[indices[1]],
            this->amplitude[indices[2]]
        );

        source->indices = indices;
        source->wavelength_index = indices[0];
        return source;
    }

    std::shared_ptr<BaseSource> get_source_by_index_sequential(const size_t index) const override {

        std::shared_ptr<BaseSource> source = std::make_shared<Planewave>(
            this->wavelength[index],
            this->jones_vector[index],
            this->amplitude[index]
        );

        return source;

    }

    void validate_sequential_data(const size_t expected_size) const override {
        // Check each vector's size and throw an error with the specific vector name if sizes don't match
        this->check_size(this->wavelength, expected_size, "wavelength");
        this->check_size(this->jones_vector, expected_size, "jones_vector");
        this->check_size(this->amplitude, expected_size, "amplitude");
    }
};
