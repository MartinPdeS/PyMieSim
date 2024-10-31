
#pragma once

#include <vector>
#include <complex>

typedef std::complex<double> complex128;

class ScattererProperties {
    private:
        std::optional<std::vector<complex128>> index_properties;
        std::optional<std::vector<std::vector<complex128>>> material_properties;

    public:
        ScattererProperties() = default;

        ScattererProperties(const std::vector<complex128> index_properties) : index_properties(index_properties) {}
        ScattererProperties(const std::vector<std::vector<complex128>> material_properties) : material_properties(material_properties) {}

        size_t size() const {
            if (index_properties)
                return index_properties->size();

            if (material_properties)
                return material_properties->size();

            throw std::runtime_error("Object not properly initialized with a valid vector.");
        }

        complex128 get (const size_t &index, const size_t &wl_index) const {
            if (index_properties)
                return (*index_properties)[index];

            if (material_properties)
                return (*material_properties)[index][wl_index];

            throw std::runtime_error("Object not properly initialized with a valid vector.");
        }
};



class MediumProperties {
    private:
        std::optional<std::vector<double>> index_properties;
        std::optional<std::vector<std::vector<double>>> material_properties;

    public:
        MediumProperties() = default;

        MediumProperties(const std::vector<double> index_properties) : index_properties(index_properties) {}
        MediumProperties(const std::vector<std::vector<double>> material_properties) : material_properties(material_properties) {}

        size_t size() const {
            if (index_properties)
                return index_properties->size();

            if (material_properties)
                return material_properties->size();

            throw std::runtime_error("Object not properly initialized with a valid vector.");
        }

        double get (const int &index, const size_t &wl_index) const {
            if (index_properties)
                return (*index_properties)[index];

            if (material_properties)
                return (*material_properties)[index][wl_index];

            throw std::runtime_error("Object not properly initialized with a valid vector.");
        }
};