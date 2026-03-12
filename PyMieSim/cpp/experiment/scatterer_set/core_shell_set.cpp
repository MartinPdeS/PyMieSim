#include "./core_shell_set.h"



// CoreShellSet methods
void CoreShellSet::update_shape() {
    this->shape = {
        core_diameter.size(),
        shell_thickness.size(),
        core_material.size(),
        shell_material.size(),
        medium.size()
    };

    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

CoreShell CoreShellSet::get_scatterer_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = this->calculate_indices(flat_index);

    CoreShell scatterer(
        core_diameter[indices[0]],
        shell_thickness[indices[1]],
        core_material[indices[2]],
        shell_material[indices[3]],
        medium[indices[4]]
    );

    scatterer.indices = indices;

    return scatterer;
}

std::shared_ptr<BaseScatterer> CoreShellSet::get_scatterer_ptr_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = this->calculate_indices(flat_index);

    CoreShell scatterer(
        core_diameter[indices[0]],
        shell_thickness[indices[1]],
        core_material[indices[2]],
        shell_material[indices[3]],
        medium[indices[4]]
    );

    scatterer.indices = indices;

    return std::make_unique<CoreShell>(scatterer);
}

void CoreShellSet::validate_sequential_data(const size_t expected_size) const {
    if (this->core_diameter.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. core_diameter has a different size than expected size.");

    if (this->shell_thickness.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. shell_thickness has a different size than expected size.");

    if (this->core_material.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. core_material has a different size than expected size.");

    if (this->shell_material.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. shell_material has a different size than expected size.");

    if (this->medium.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. medium has a different size than expected size.");
}

CoreShell CoreShellSet::get_scatterer_by_index_sequential(const size_t index) const {
    CoreShell scatterer(
        core_diameter[index],
        shell_thickness[index],
        core_material[index],
        shell_material[index],
        medium[index]
    );

    return scatterer;
}

std::shared_ptr<BaseScatterer> CoreShellSet::get_scatterer_ptr_by_index_sequential(const size_t index) const {
    CoreShell scatterer = CoreShell(
        core_diameter[index],
        shell_thickness[index],
        core_material[index],
        shell_material[index],
        medium[index]
    );

    return std::make_unique<CoreShell>(scatterer);
}
