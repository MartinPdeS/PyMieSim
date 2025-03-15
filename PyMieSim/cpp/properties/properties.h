
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

        ScattererProperties(const std::vector<complex128> index_properties);
        ScattererProperties(const std::vector<std::vector<complex128>> material_properties);

        size_t size() const;

        complex128 get (const size_t &index, const size_t &wl_index) const;
};



class MediumProperties {
    private:
        std::optional<std::vector<double>> index_properties;
        std::optional<std::vector<std::vector<double>>> material_properties;

    public:
        MediumProperties() = default;

        MediumProperties(const std::vector<double> index_properties);
        MediumProperties(const std::vector<std::vector<double>> material_properties);

        size_t size() const;

        double get (const int &index, const size_t &wl_index) const;
};