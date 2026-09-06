#include "./core_shell_set.h"

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

std::shared_ptr<BaseScatterer> CoreShellSet::get_scatterer_by_index(const size_t flat_index) const {
    std::vector<size_t> indices = this->calculate_indices(flat_index);

    std::shared_ptr<CoreShell> scatterer = std::make_shared<CoreShell>(
        core_diameter[indices[0]],
        shell_thickness[indices[1]],
        core_material[indices[2]]->clone(),
        shell_material[indices[3]]->clone(),
        medium[indices[4]]->clone()
    );

    scatterer->indices = indices;

    return scatterer;
}

std::shared_ptr<BaseScatterer> CoreShellSet::get_scatterer_by_index_sequential(const size_t index) const {

    std::shared_ptr<CoreShell> scatterer = std::make_shared<CoreShell>(
        core_diameter[index],
        shell_thickness[index],
        core_material[index]->clone(),
        shell_material[index]->clone(),
        medium[index]->clone()
    );

    return scatterer;
}
