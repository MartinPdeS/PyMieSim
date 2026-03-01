#include "./source_set.h"


std::shared_ptr<BaseSource>
GaussianSourceSet::get_source_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    std::shared_ptr<BaseSource> source = std::make_shared<Gaussian>(
        this->wavelength[indices[0]],
        this->polarization->polarization_states[indices[1]],
        this->numerical_aperture[indices[2]],
        this->optical_power[indices[3]]
    );

    source->indices = indices;
    source->wavelength_index = indices[0];

    return source;
}

std::shared_ptr<BaseSource> GaussianSourceSet::get_source_by_index_sequential(const size_t index) const {
    std::shared_ptr<BaseSource> source = std::make_shared<Gaussian>(
        this->wavelength[index],
        this->polarization->polarization_states[index],
        this->numerical_aperture[index],
        this->optical_power[index]
    );

    return source;
}

void GaussianSourceSet::validate_sequential_data(const size_t expected_size) const {
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    this->check_size(this->wavelength, expected_size, "wavelength");
    this->check_size(this->polarization->polarization_states, expected_size, "polarization");
    this->check_size(this->numerical_aperture, expected_size, "numerical_aperture");
    this->check_size(this->optical_power, expected_size, "optical_power");
}


std::shared_ptr<BaseSource> PlaneWaveSourceSet::get_source_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    std::shared_ptr<BaseSource> source = std::make_shared<Planewave>(
        this->wavelength[indices[0]],
        this->polarization->polarization_states[indices[1]],
        this->amplitude[indices[2]]
    );

    source->indices = indices;
    source->wavelength_index = indices[0];
    return source;
}

std::shared_ptr<BaseSource> PlaneWaveSourceSet::get_source_by_index_sequential(const size_t index) const {

    std::shared_ptr<BaseSource> source = std::make_shared<Planewave>(
        this->wavelength[index],
        this->polarization->polarization_states[index],
        this->amplitude[index]
    );

    return source;

}

void PlaneWaveSourceSet::validate_sequential_data(const size_t expected_size) const {
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    this->check_size(this->wavelength, expected_size, "wavelength");
    this->check_size(this->polarization->polarization_states, expected_size, "polarization");
    this->check_size(this->amplitude, expected_size, "amplitude");
}

