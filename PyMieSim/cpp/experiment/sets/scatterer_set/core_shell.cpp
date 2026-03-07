#include "./core_shell.h"



// CoreShellSet methods
void CoreShellSet::update_shape() {
    this->shape = {
        core_diameter.size(),
        shell_thickness.size(),
        core_property.size(),
        shell_property.size(),
        medium_property.size()
    };

    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

CoreShell CoreShellSet::get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = this->calculate_indices(flat_index);

    CoreShell scatterer(
        core_diameter[indices[0]],
        shell_thickness[indices[1]],
        core_property.get(indices[2], source->wavelength_index),
        shell_property.get(indices[3], source->wavelength_index),
        medium_property.get(indices[4], source->wavelength_index),
        source
    );

    scatterer.indices = indices;

    return scatterer;
}

std::unique_ptr<BaseScatterer> CoreShellSet::get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = this->calculate_indices(flat_index);

    CoreShell scatterer(
        core_diameter[indices[0]],
        shell_thickness[indices[1]],
        core_property.get(indices[2], source->wavelength_index),
        shell_property.get(indices[3], source->wavelength_index),
        medium_property.get(indices[4], source->wavelength_index),
        source
    );

    scatterer.indices = indices;

    return std::make_unique<CoreShell>(scatterer);
}

void CoreShellSet::validate_sequential_data(const size_t expected_size) const {
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    if (this->core_diameter.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. core_diameter has a different size than expected size.");

    if (this->shell_thickness.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. shell_thickness has a different size than expected size.");

    if (this->core_property.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. core_property has a different size than expected size.");

    if (this->shell_property.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. shell_property has a different size than expected size.");

    if (this->medium_property.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. medium_property has a different size than expected size.");
}

CoreShell CoreShellSet::get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    CoreShell scatterer(
        core_diameter[index],
        shell_thickness[index],
        core_property.get(index, source->wavelength_index),
        shell_property.get(index, source->wavelength_index),
        medium_property.get(index, source->wavelength_index),
        source
    );

    return scatterer;
}

std::unique_ptr<BaseScatterer> CoreShellSet::get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    CoreShell scatterer = CoreShell(
        core_diameter[index],
        shell_thickness[index],
        core_property.get(index, source->wavelength_index),
        shell_property.get(index, source->wavelength_index),
        medium_property.get(index, source->wavelength_index),
        source
    );

    return std::make_unique<CoreShell>(scatterer);
}
