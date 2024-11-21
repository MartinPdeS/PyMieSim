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
        MediumProperties medium_property;

        Set() = default;

        Set(const std::vector<double>& diameter, const ScattererProperties& property, const MediumProperties& medium_property)
            : diameter(diameter), property(property), medium_property(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            this->shape = {
                diameter.size(),
                property.size(),
                medium_property.size()
            };
        }

        void validate_sequential_data(const size_t expected_size) const {
            // Check each vector's size and throw an error with the specific vector name if sizes don't match
            if (this->diameter.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. diameter has a different size than expected size.");

            if (this->property.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. property has a different size than expected size.");

            if (this->medium_property.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. medium_property has a different size than expected size.");
        }

        Scatterer get_scatterer_by_index_sequential(const size_t index, const SOURCE::BaseSource& source) const {
            return Scatterer(
                this->diameter[index],
                this->property.get(index, source.wavelength_index),
                this->medium_property.get(index, source.wavelength_index),
                source
            );
        }

        Scatterer get_scatterer_by_index(const size_t flat_index, const SOURCE::BaseSource& source) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            Scatterer scatterer(
                this->diameter[indices[0]],
                this->property.get(indices[1], source.wavelength_index),
                this->medium_property.get(indices[2], source.wavelength_index),
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
        MediumProperties medium_property;


        Set() = default;
        Set(const std::vector<double>& diameter, const ScattererProperties& property, const MediumProperties& medium_property)
            : diameter(diameter), property(property), medium_property(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            this->shape = {
                diameter.size(),
                property.size(),
                medium_property.size()
            };
        }

        void validate_sequential_data(const size_t expected_size) const {
            // Check each vector's size and throw an error with the specific vector name if sizes don't match
            if (this->diameter.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. diameter has a different size than expected size.");

            if (this->property.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. property has a different size than expected size.");

            if (this->medium_property.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. medium_property has a different size than expected size.");
        }

        Scatterer get_scatterer_by_index_sequential(const size_t index, const SOURCE::BaseSource& source) const {
            return Scatterer(
                this->diameter[index],
                this->property.get(index, source.wavelength_index),
                this->medium_property.get(index, source.wavelength_index),
                source
            );
        }

        Scatterer get_scatterer_by_index(const size_t flat_index, const SOURCE::BaseSource& source) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            Scatterer scatterer(
                diameter[indices[0]],
                property.get(indices[1], source.wavelength_index),
                medium_property.get(indices[2], source.wavelength_index),
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
        std::vector<double> shell_thickness;
        ScattererProperties core_property;
        ScattererProperties shell_property;
        MediumProperties medium_property;

        Set() = default;

        Set(
            const std::vector<double>& core_diameter,
            const std::vector<double>& shell_thickness,
            const ScattererProperties& core_property,
            const ScattererProperties& shell_property,
            const MediumProperties& medium_property)
            : core_diameter(core_diameter), shell_thickness(shell_thickness), core_property(core_property), shell_property(shell_property), medium_property(medium_property)
            {
                update_shape();
                total_combinations = get_vector_sigma(shape);
            }

        void update_shape() override {
            this->shape = {
                core_diameter.size(),
                shell_thickness.size(),
                core_property.size(),
                shell_property.size(),
                medium_property.size()
            };
        }

        Scatterer get_scatterer_by_index(const size_t flat_index, SOURCE::BaseSource& source) const {

            std::vector<size_t> indices = this->calculate_indices(flat_index);

            Scatterer scatterer(
                core_diameter[indices[0]],
                shell_thickness[indices[1]],
                core_property.get(indices[2], source.wavelength_index),
                shell_property.get(indices[3], source.wavelength_index),
                medium_property.get(indices[4], source.wavelength_index),
                source
            );

            scatterer.indices = indices;

            return scatterer;
        }

        void validate_sequential_data(const size_t expected_size) const {
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

        Scatterer get_scatterer_by_index_sequential(const size_t index, SOURCE::BaseSource& source) const {
            Scatterer scatterer(
                core_diameter[index],
                shell_thickness[index],
                core_property.get(index, source.wavelength_index),
                shell_property.get(index, source.wavelength_index),
                medium_property.get(index, source.wavelength_index),
                source
            );

            return scatterer;
        }
    };
}