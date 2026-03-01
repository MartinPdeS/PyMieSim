#include "./scatterer_set.h"



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
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    this->check_size(this->diameter, expected_size, "diameter");
    this->check_size(this->property, expected_size, "property");
    this->check_size(this->medium_property, expected_size, "medium_property");
}

Sphere SphereSet::get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    return Sphere(
        this->diameter[index],
        this->property.get(index, source->wavelength_index),
        this->medium_property.get(index, source->wavelength_index),
        source
    );
}

std::unique_ptr<BaseScatterer>
SphereSet::get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    Sphere scatterer(
        this->diameter[index],
        this->property.get(index, source->wavelength_index),
        this->medium_property.get(index, source->wavelength_index),
        source
    );

    return std::make_unique<Sphere>(scatterer);
}

Sphere SphereSet::get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    Sphere scatterer(
        this->diameter[indices[0]],
        this->property.get(indices[1], source->wavelength_index),
        this->medium_property.get(indices[2], source->wavelength_index),
        source
    );

    scatterer.indices = indices;

    return scatterer;
}

std::unique_ptr<BaseScatterer> SphereSet::get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    Sphere scatterer(
        this->diameter[indices[0]],
        this->property.get(indices[1], source->wavelength_index),
        this->medium_property.get(indices[2], source->wavelength_index),
        source
    );

    scatterer.indices = indices;

    return std::make_unique<Sphere>(scatterer);
}




// InfiniteCylinderSet methods
void InfiniteCylinderSet::update_shape() {
    this->shape = {
        diameter.size(),
        property.size(),
        medium_property.size()
    };

    total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
}

void InfiniteCylinderSet::validate_sequential_data(const size_t expected_size) const {
    // Check each vector's size and throw an error with the specific vector name if sizes don't match
    if (this->diameter.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. diameter has a different size than expected size.");

    if (this->property.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. property has a different size than expected size.");

    if (this->medium_property.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. medium_property has a different size than expected size.");
}

InfiniteCylinder InfiniteCylinderSet::get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    return InfiniteCylinder(
        this->diameter[index],
        this->property.get(index, source->wavelength_index),
        this->medium_property.get(index, source->wavelength_index),
        source
    );
}

std::unique_ptr<BaseScatterer> InfiniteCylinderSet::get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
    InfiniteCylinder scatterer(
        this->diameter[index],
        this->property.get(index, source->wavelength_index),
        this->medium_property.get(index, source->wavelength_index),
        source
    );

    return std::make_unique<InfiniteCylinder>(scatterer);
}

InfiniteCylinder InfiniteCylinderSet::get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    InfiniteCylinder scatterer(
        diameter[indices[0]],
        property.get(indices[1], source->wavelength_index),
        medium_property.get(indices[2], source->wavelength_index),
        source
    );

    scatterer.indices = indices;

    return scatterer;
}


std::unique_ptr<BaseScatterer> InfiniteCylinderSet::get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
    std::vector<size_t> indices = calculate_indices(flat_index);

    InfiniteCylinder scatterer = InfiniteCylinder(
        diameter[indices[0]],
        property.get(indices[1], source->wavelength_index),
        medium_property.get(indices[2], source->wavelength_index),
        source
    );

    scatterer.indices = indices;

    return std::make_unique<InfiniteCylinder>(scatterer);
}

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
