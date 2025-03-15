#include "sets/base.h"


BaseSet::BaseSet(bool is_sequential){this->is_sequential = is_sequential;}

// Calculate the multi-dimensional indices for the current index
std::vector<size_t>
BaseSet::calculate_indices(size_t flat_index) const
{
    std::vector<size_t> indices(shape.size());
    for (size_t i = shape.size(); i-- > 0;)
    {
        indices[i] = flat_index % shape[i];
        flat_index /= shape[i];
    }

    return indices;
}

template <typename Container>
void BaseSet::check_size(const Container& container, size_t expected_size, const std::string& name) const {
    if (container.size() != expected_size)
        throw std::runtime_error("Error: Vector size mismatch in sequential computation. " + name + " has size " + std::to_string(container.size()) + ", but wavelength size is " + std::to_string(expected_size) + ".");
}

