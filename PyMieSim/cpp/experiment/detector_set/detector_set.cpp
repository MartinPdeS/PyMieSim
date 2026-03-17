#include "./detector_set.h"


PhotodiodeSet::PhotodiodeSet(
    const std::vector<unsigned> &sampling,
    const std::vector<double> &numerical_aperture,
    const std::vector<double> &cache_numerical_aperture,
    const std::vector<double> &phi_offset,
    const std::vector<double> &gamma_offset,
    const PolarizationSet &polarization_filter_set,
    const MediumSet &medium_set,
    const bool is_sequential
)
:   sampling(sampling),
    numerical_aperture(numerical_aperture),
    cache_numerical_aperture(cache_numerical_aperture),
    phi_offset(phi_offset),
    gamma_offset(gamma_offset),
    polarization_filter_set(polarization_filter_set),
    medium(medium_set)
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
        polarization_filter_set.number_of_states(),
        medium.size()
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
        this->polarization_filter_set.polarization_states[indices[5]],
        this->medium[indices[6]]
    );

    detector->indices = indices;

    return detector;
}

std::shared_ptr<BaseDetector>
PhotodiodeSet::get_detector_by_index_sequential(size_t index) const {
    std::shared_ptr<Photodiode> detector = std::make_shared<Photodiode>(
        this->sampling[index],
        this->numerical_aperture[index],
        this->cache_numerical_aperture[index],
        this->phi_offset[index],
        this->gamma_offset[index],
        this->polarization_filter_set.polarization_states[index],
        this->medium[index]
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
    const PolarizationSet &polarization_filter_set,
    const std::vector<double> &rotation,
    const MediumSet &medium_set,
    const bool &mean_coupling,
    const bool is_sequential
)
:   mode_numbers(mode_numbers),
    sampling(sampling),
    numerical_aperture(numerical_aperture),
    cache_numerical_aperture(cache_numerical_aperture),
    phi_offset(phi_offset),
    gamma_offset(gamma_offset),
    polarization_filter_set(polarization_filter_set),
    rotation(rotation),
    medium(medium_set),
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
        polarization_filter_set.number_of_states(),
        rotation.size(),
        medium.size()
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
        this->polarization_filter_set.polarization_states[indices[6]],
        this->rotation[indices[7]],
        this->mean_coupling,
        this->medium[indices[8]]
    );

    detector->indices = indices;

    return detector;
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
        this->polarization_filter_set.polarization_states[index],
        this->rotation[index],
        this->mean_coupling,
        this->medium[index]
    );

    detector->indices = {index};
    return detector;
}
