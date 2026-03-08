#include "./detector_set.h"


PhotodiodeSet::PhotodiodeSet(
    const std::vector<unsigned> &sampling,
    const std::vector<double> &numerical_aperture,
    const std::vector<double> &cache_numerical_aperture,
    const std::vector<double> &phi_offset,
    const std::vector<double> &gamma_offset,
    const std::vector<double> &polarization_filter,
    const std::vector<double> &medium_refractive_index,
    const bool is_sequential
)
:   sampling(sampling),
    numerical_aperture(numerical_aperture),
    cache_numerical_aperture(cache_numerical_aperture),
    phi_offset(phi_offset),
    gamma_offset(gamma_offset),
    polarization_filter(polarization_filter),
    medium_refractive_index(medium_refractive_index)
    {
        this->is_sequential = is_sequential;
        this->is_empty = false;
        this->update_shape();
    }

void PhotodiodeSet::update_shape() {
    this->shape = {
        sampling.size(),
        numerical_aperture.size(),
        cache_numerical_aperture.size(),
        phi_offset.size(),
        gamma_offset.size(),
        polarization_filter.size(),
        medium_refractive_index.size()
    };
    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

std::shared_ptr<BaseDetector>
PhotodiodeSet::get_detector_by_index(long long flat_index) const {
    std::vector<size_t> indices = this->calculate_indices(flat_index);

    std::shared_ptr<Photodiode> detector = std::make_shared<Photodiode>(
        this->sampling[indices[0]],
        this->numerical_aperture[indices[1]],
        this->cache_numerical_aperture[indices[2]],
        this->phi_offset[indices[3]],
        this->gamma_offset[indices[4]],
        this->polarization_filter[indices[5]],
        this->medium_refractive_index[indices[6]]
    );

    detector->indices = indices;

    return detector;
}

void PhotodiodeSet::validate_sequential_data(const size_t expected_size) const {
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    this->check_size(this->sampling, expected_size, "sampling");
    this->check_size(this->numerical_aperture, expected_size, "numerical_aperture");
    this->check_size(this->cache_numerical_aperture, expected_size, "cache_numerical_aperture");

    this->check_size(this->phi_offset, expected_size, "phi_offset");
    this->check_size(this->gamma_offset, expected_size, "gamma_offset");
    this->check_size(this->polarization_filter, expected_size, "polarization_filter");
    this->check_size(this->medium_refractive_index, expected_size, "medium_refractive_index");
}

std::shared_ptr<BaseDetector>
PhotodiodeSet::get_detector_by_index_sequential(size_t index) const {
    std::shared_ptr<Photodiode> detector = std::make_shared<Photodiode>(
        this->sampling[index],
        this->numerical_aperture[index],
        this->cache_numerical_aperture[index],
        this->phi_offset[index],
        this->gamma_offset[index],
        this->polarization_filter[index],
        this->medium_refractive_index[index]
    );
    detector->indices = {index};
    return detector;
}





CoherentModeSet::CoherentModeSet(
    const std::vector<std::string> &mode_numbers,
    const std::vector<unsigned> &sampling,
    const std::vector<double> &numerical_aperture,
    const std::vector<double> &cache_numerical_aperture,
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
    numerical_aperture(numerical_aperture),
    cache_numerical_aperture(cache_numerical_aperture),
    phi_offset(phi_offset),
    gamma_offset(gamma_offset),
    polarization_filter(polarization_filter),
    rotation(rotation),
    medium_refractive_index(medium_refractive_index),
    mean_coupling(mean_coupling)
    {
        for (std::string mode_number : mode_numbers) {
            if (mode_number.size() < 4 || !std::isalpha(mode_number[0]) || !std::isalpha(mode_number[1]) || !std::isdigit(mode_number[2]) || !std::isdigit(mode_number[3])) {
                throw std::invalid_argument("Invalid mode number format: " + mode_number + ". Expected format like 'LP01', 'HG12', etc.");
            }
        }
        this->is_sequential = is_sequential;
        this->is_empty = false;
        this->update_shape();
    }

void CoherentModeSet::update_shape() {
    this->shape = {
        mode_numbers.size(),
        sampling.size(),
        numerical_aperture.size(),
        cache_numerical_aperture.size(),
        phi_offset.size(),
        gamma_offset.size(),
        polarization_filter.size(),
        rotation.size(),
        medium_refractive_index.size()
    };
    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

std::shared_ptr<BaseDetector>
CoherentModeSet::get_detector_by_index(long long flat_index) const {
    std::vector<size_t> indices = this->calculate_indices(flat_index);
    std::shared_ptr<CoherentMode> detector = std::make_shared<CoherentMode>(
        this->mode_numbers[indices[0]],
        this->sampling[indices[1]],
        this->numerical_aperture[indices[2]],
        this->cache_numerical_aperture[indices[3]],
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

void CoherentModeSet::validate_sequential_data(const size_t expected_size) const {
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    this->check_size(this->mode_numbers, expected_size, "mode_numbers");
    this->check_size(this->sampling, expected_size, "sampling");
    this->check_size(this->numerical_aperture, expected_size, "numerical_aperture");
    this->check_size(this->cache_numerical_aperture, expected_size, "cache_numerical_aperture");

    this->check_size(this->phi_offset, expected_size, "phi_offset");
    this->check_size(this->gamma_offset, expected_size, "gamma_offset");
    this->check_size(this->polarization_filter, expected_size, "polarization_filter");
    this->check_size(this->rotation, expected_size, "rotation");
    this->check_size(this->medium_refractive_index, expected_size, "medium_refractive_index");
}

std::shared_ptr<BaseDetector>
CoherentModeSet::get_detector_by_index_sequential(size_t index) const {
    std::shared_ptr<CoherentMode> detector = std::make_shared<CoherentMode>(
        this->mode_numbers[index],
        this->sampling[index],
        this->numerical_aperture[index],
        this->cache_numerical_aperture[index],
        this->phi_offset[index],
        this->gamma_offset[index],
        this->polarization_filter[index],
        this->rotation[index],
        this->mean_coupling,
        this->medium_refractive_index[index]
    );

    detector->indices = {index};
    return detector;
}
