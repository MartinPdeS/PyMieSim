#include "./sphere_set.h"



// SphereSet methods
void SphereSet::update_shape() {
    this->shape = {
        diameter.size(),
        material.size(),
        medium.size()
    };

    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

void SphereSet::validate_sequential_data(const size_t expected_size) const {
    this->check_size(this->diameter, expected_size, "diameter");
    this->check_size(this->material, expected_size, "material");
    this->check_size(this->medium, expected_size, "medium");
}

std::shared_ptr<BaseScatterer>
SphereSet::get_scatterer_by_index_sequential(const size_t index) const {

    std::shared_ptr<Sphere> scatterer = std::make_shared<Sphere>(
        this->diameter[index],
        this->material[index],
        this->medium[index]
    );

    return scatterer;


}

std::shared_ptr<BaseScatterer>
SphereSet::get_scatterer_ptr_by_index_sequential(const size_t index) const {

    std::shared_ptr<Sphere> scatterer = std::make_shared<Sphere>(
        this->diameter[index],
        this->material[index],
        this->medium[index]
    );

    return scatterer;
}

std::shared_ptr<BaseScatterer>
SphereSet::get_scatterer_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    std::shared_ptr<Sphere> scatterer = std::make_shared<Sphere>(
        this->diameter[indices[0]],
        this->material[indices[1]],
        this->medium[indices[2]]
    );

    scatterer->indices = indices;

    return scatterer;
}

std::shared_ptr<BaseScatterer> SphereSet::get_scatterer_ptr_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    std::shared_ptr<Sphere> scatterer = std::make_shared<Sphere>(
        this->diameter[indices[0]],
        this->material[indices[1]],
        this->medium[indices[2]]
    );

    scatterer->indices = indices;

    return scatterer;
}
