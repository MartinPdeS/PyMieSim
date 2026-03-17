#include "./source_set.h"


std::shared_ptr<BaseSource>
GaussianSourceSet::get_source_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    std::shared_ptr<BaseSource> source = std::make_shared<Gaussian>(
        this->wavelength[indices[0]],
        this->polarization.polarization_states[indices[1]],
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
        this->polarization.polarization_states[index],
        this->numerical_aperture[index],
        this->optical_power[index]
    );

    return source;
}


std::shared_ptr<BaseSource> PlaneWaveSourceSet::get_source_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    std::shared_ptr<BaseSource> source = std::make_shared<Planewave>(
        this->wavelength[indices[0]],
        this->polarization.polarization_states[indices[1]],
        this->amplitude[indices[2]]
    );

    source->indices = indices;
    source->wavelength_index = indices[0];
    return source;
}

std::shared_ptr<BaseSource> PlaneWaveSourceSet::get_source_by_index_sequential(const size_t index) const {

    std::shared_ptr<BaseSource> source = std::make_shared<Planewave>(
        this->wavelength[index],
        this->polarization.polarization_states[index],
        this->amplitude[index]
    );

    return source;

}
