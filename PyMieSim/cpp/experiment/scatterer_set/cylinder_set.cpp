#include "./cylinder_set.h"


// InfiniteCylinderSet methods
void InfiniteCylinderSet::update_shape() {
    this->shape = {
        diameter.size(),
        material.size(),
        medium.size()
    };

    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

void InfiniteCylinderSet::validate_sequential_data(const size_t expected_size) const {
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    if (this->diameter.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. diameter has a different size than expected size.");

    if (this->material.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. material has a different size than expected size.");

    if (this->medium.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. medium has a different size than expected size.");
}

std::shared_ptr<BaseScatterer> InfiniteCylinderSet::get_scatterer_by_index_sequential(const size_t index) const {

    std::shared_ptr<InfiniteCylinder> scatterer = std::make_shared<InfiniteCylinder>(
        this->diameter[index],
        this->material[index]->clone(),
        this->medium[index]->clone()
    );

    return scatterer;
}


std::shared_ptr<BaseScatterer> InfiniteCylinderSet::get_scatterer_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    std::shared_ptr<InfiniteCylinder> scatterer = std::make_shared<InfiniteCylinder>(
        diameter[indices[0]],
        material[indices[1]]->clone(),
        medium[indices[2]]->clone()
    );

    scatterer->indices = indices;

    return scatterer;
}
