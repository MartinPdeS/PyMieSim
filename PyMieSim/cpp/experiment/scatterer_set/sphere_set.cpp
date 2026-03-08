#include "./sphere_set.h"



// SphereSet methods
void SphereSet::update_shape() {
    this->shape = {
        diameter.size(),
        property.size(),
        medium_property.size()
    };

    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

void SphereSet::validate_sequential_data(const size_t expected_size) const {
    this->check_size(this->diameter, expected_size, "diameter");
    this->check_size(this->property, expected_size, "property");
    this->check_size(this->medium_property, expected_size, "medium_property");
}

Sphere SphereSet::get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    return Sphere(
        this->diameter[index],
        this->property[index],
        this->medium_property[index],
        source
    );
}

std::unique_ptr<BaseScatterer>
SphereSet::get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    Sphere scatterer(
        this->diameter[index],
        this->property[index],
        this->medium_property[index],
        source
    );

    return std::make_unique<Sphere>(scatterer);
}

Sphere SphereSet::get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    Sphere scatterer(
        this->diameter[indices[0]],
        this->property[indices[1]],
        this->medium_property[indices[2]],
        source
    );

    scatterer.indices = indices;

    return scatterer;
}

std::unique_ptr<BaseScatterer> SphereSet::get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    Sphere scatterer(
        this->diameter[indices[0]],
        this->property[indices[1]],
        this->medium_property[indices[2]],
        source
    );

    scatterer.indices = indices;

    return std::make_unique<Sphere>(scatterer);
}
