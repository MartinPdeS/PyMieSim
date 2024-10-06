#pragma once

#include "experiment/includes/sets/base.cpp"
#include "single/includes/sphere.cpp"
#include "single/includes/cylinder.cpp"
#include "single/includes/coreshell.cpp"

using complex128 = std::complex<double>;

// Sphere class inheriting from BaseSet
namespace SPHERE {
    class Set : public BaseSet
    {
    public:
        std::vector<double> diameter;
        std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> property;
        std::variant<std::vector<double>, std::vector<std::vector<double>>> medium;

        Set() = default;

        Set(
            const std::vector<double>& diameter,
            const std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>>& property,
            const std::variant<std::vector<double>, std::vector<std::vector<double>>>& medium_property)
            : diameter(diameter), property(property), medium(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            shape.clear();
            shape.push_back(diameter.size());
            shape.push_back(get_variant_size(property));
            shape.push_back(get_variant_size(medium));
        }

        Scatterer get_scatterer_by_index(size_t flat_index, SOURCE::BaseSource& source) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            Scatterer scatterer(
                diameter[indices[0]],
                get_variant_value(property, indices[1], source.indices[0]),
                get_variant_value(medium, indices[2], source.indices[0]),
                source
            );

            scatterer.indices = indices;

            return scatterer;
        }
    };
}

// Cylinder class inheriting from BaseSet
namespace CYLINDER {
    class Set : public BaseSet {
    public:
        std::vector<double> diameter;
        std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> property;
        std::variant<std::vector<double>, std::vector<std::vector<double>>> medium;

        Set() = default;
        Set(
            const std::vector<double>& diameter,
            const std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>>& property,
            const std::variant<std::vector<double>, std::vector<std::vector<double>>>& medium_property)
            : diameter(diameter), property(property), medium(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override
        {
            shape.clear();
            shape.push_back(diameter.size());
            shape.push_back(get_variant_size(property));
            shape.push_back(get_variant_size(medium));
        }

        Scatterer get_scatterer_by_index(size_t flat_index, SOURCE::BaseSource& source) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            Scatterer scatterer(
                diameter[indices[0]],
                get_variant_value(property, indices[1], source.indices[0]),
                get_variant_value(medium, indices[2], source.indices[0]),
                source
            );

            scatterer.indices = indices;

            return scatterer;
        }
    };
}

// Core-shell class
namespace CORESHELL {
    class Set : public BaseSet {
    public:
        std::vector<double> core_diameter;
        std::vector<double> shell_width;
        std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> core_property;
        std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> shell_property;
        std::variant<std::vector<double>, std::vector<std::vector<double>>> medium;

        Set() = default;

        Set(
            const std::vector<double>& core_diameter,
            const std::vector<double>& shell_width,
            const std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>>& core_property,
            const std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>>& shell_property,
            const std::variant<std::vector<double>, std::vector<std::vector<double>>>& medium_property)
            : core_diameter(core_diameter), shell_width(shell_width), core_property(core_property), shell_property(shell_property), medium(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            shape.clear();
            shape.push_back(core_diameter.size());
            shape.push_back(shell_width.size());
            shape.push_back(get_variant_size(core_property));
            shape.push_back(get_variant_size(shell_property));
            shape.push_back(get_variant_size(medium));
        }

        Scatterer get_scatterer_by_index(size_t flat_index, SOURCE::BaseSource& source) const {

            std::vector<size_t> indices = this->calculate_indices(flat_index);

            Scatterer scatterer(
                core_diameter[indices[0]],
                shell_width[indices[1]],
                get_variant_value(core_property, indices[2], source.indices[0]),
                get_variant_value(shell_property, indices[3], source.indices[0]),
                get_variant_value(medium, indices[4], source.indices[0]),
                source
            );

            scatterer.indices = indices;

            return scatterer;
        }
    };
}