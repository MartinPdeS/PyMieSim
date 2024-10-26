
#pragma once

#include <vector>
#include <complex>

typedef std::complex<double> complex128;

class ScattererProperties {
    public:
        std::vector<complex128> c_index_properties;
        std::vector<std::vector<complex128>> c_material_properties;
        std::vector<double> d_index_properties;
        std::vector<std::vector<double>> d_material_properties;

        ScattererProperties(std::vector<complex128> index_properties) : c_index_properties(index_properties) {}
        ScattererProperties(std::vector<std::vector<complex128>> material_properties) : c_material_properties(material_properties) {}
        ScattererProperties(std::vector<double> index_properties) : d_index_properties(index_properties) {}
        ScattererProperties(std::vector<std::vector<double>> material_properties) : d_material_properties(material_properties) {}

        size_t size()
        {
            if (c_index_properties.size() != 0)
                return c_index_properties.size();

            if (c_material_properties.size() != 0)
                return c_material_properties.size();

            if (d_index_properties.size() != 0)
                return d_index_properties.size();

            if (d_material_properties.size() != 0)
                return d_material_properties.size();
        }
};