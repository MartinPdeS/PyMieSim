#pragma once

#include "experiment/includes/scatterer_properties.cpp"
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
        ScattererProperties property;
        MediumProperties medium;

        Set() = default;

        Set(
            const std::vector<double>& diameter,
            const ScattererProperties& property,
            const MediumProperties& medium_property)
            : diameter(diameter), property(property), medium(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            this->shape = {
                diameter.size(),
                property.size(),
                medium.size()
            };
        }

        Scatterer get_scatterer_by_index(size_t flat_index, SOURCE::BaseSource& source) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            Scatterer scatterer(
                diameter[indices[0]],
                property.get(indices[1], source.indices[0]),
                medium.get(indices[2], source.indices[0]),
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
        ScattererProperties property;
        MediumProperties medium;


        Set() = default;
        Set(
            const std::vector<double>& diameter,
            const ScattererProperties& property,
            const MediumProperties& medium_property)
            : diameter(diameter), property(property), medium(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            this->shape = {
                diameter.size(),
                property.size(),
                medium.size()
            };
        }

        Scatterer get_scatterer_by_index(size_t flat_index, SOURCE::BaseSource& source) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            Scatterer scatterer(
                diameter[indices[0]],
                property.get(indices[1], source.indices[0]),
                medium.get(indices[2], source.indices[0]),
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
        ScattererProperties core_property;
        ScattererProperties shell_property;
        MediumProperties medium;

        Set() = default;

        Set(
            const std::vector<double>& core_diameter,
            const std::vector<double>& shell_width,
            const ScattererProperties& core_property,
            const ScattererProperties& shell_property,
            const MediumProperties& medium_property)
            : core_diameter(core_diameter), shell_width(shell_width), core_property(core_property), shell_property(shell_property), medium(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            this->shape = {
                core_diameter.size(),
                shell_width.size(),
                core_property.size(),
                shell_property.size(),
                medium.size()
            };
        }

        Scatterer get_scatterer_by_index(size_t flat_index, SOURCE::BaseSource& source) const {

            std::vector<size_t> indices = this->calculate_indices(flat_index);

            Scatterer scatterer(
                core_diameter[indices[0]],
                shell_width[indices[1]],
                core_property.get(indices[2], source.indices[0]),
                shell_property.get(indices[3], source.indices[0]),
                medium.get(indices[4], source.indices[0]),
                source
            );

            scatterer.indices = indices;

            return scatterer;
        }
    };
}