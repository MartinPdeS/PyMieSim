#include <vector>
#include <complex>

#include "properties/properties.h"

typedef std::complex<double> complex128;


ScattererProperties::ScattererProperties(const std::vector<complex128> index_properties) : index_properties(index_properties) {}
ScattererProperties::ScattererProperties(const std::vector<std::vector<complex128>> material_properties) : material_properties(material_properties) {}

size_t ScattererProperties::size() const {
    if (index_properties)
        return index_properties->size();

    if (material_properties)
        return material_properties->size();

    throw std::runtime_error("Object not properly initialized with a valid vector.");
}

complex128 ScattererProperties::get(const size_t &index, const size_t &wl_index) const {
    if (index_properties)
        return (*index_properties)[index];

    if (material_properties)
        return (*material_properties)[index][wl_index];

    throw std::runtime_error("Object not properly initialized with a valid vector.");
}


MediumProperties::MediumProperties(const std::vector<double> index_properties) : index_properties(index_properties) {}
MediumProperties::MediumProperties(const std::vector<std::vector<double>> material_properties) : material_properties(material_properties) {}

size_t MediumProperties::size() const {
    if (index_properties)
        return index_properties->size();

    if (material_properties)
        return material_properties->size();

    throw std::runtime_error("Object not properly initialized with a valid vector.");
}

double MediumProperties::get(const int &index, const size_t &wl_index) const {
    if (index_properties)
        return (*index_properties)[index];

    if (material_properties)
        return (*material_properties)[index][wl_index];

    throw std::runtime_error("Object not properly initialized with a valid vector.");
}
