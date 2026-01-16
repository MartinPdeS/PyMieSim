#pragma once

#include <memory>
#include "properties.h"
#include "base_set.h"
#include "single/scatterer/base_scatterer/base_scatterer.h"

using ScattererPtr = std::unique_ptr<BaseScatterer>;   // const because callers must not mutate



class ScattererSet: public BaseSet
{
public:
    ScattererSet() = default;
    ScattererSet(const bool is_sequential) : BaseSet(is_sequential) {}

    virtual void validate_sequential_data(const size_t expected_size) const = 0;
    virtual std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const = 0;
    virtual std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const = 0;

};





// Sphere class inheriting from BaseSet
class SphereSet : public ScattererSet
{
public:
    std::vector<double> diameter;
    ScattererProperties property;
    MediumProperties medium_property;

    SphereSet() = default;

    SphereSet(const std::vector<double>& diameter, const ScattererProperties& property, const MediumProperties& medium_property, const bool is_sequential)
        : ScattererSet(is_sequential), diameter(diameter), property(property), medium_property(medium_property)
        {this->update_shape();}

    void update_shape() override {
        this->shape = {
            diameter.size(),
            property.size(),
            medium_property.size()
        };

        total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
    }

    void validate_sequential_data(const size_t expected_size) const override {
        // Check each vector's size and throw an error with the specific vector name if sizes don't match
        this->check_size(this->diameter, expected_size, "diameter");
        this->check_size(this->property, expected_size, "property");
        this->check_size(this->medium_property, expected_size, "medium_property");
    }

    Sphere get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
        return Sphere(
            this->diameter[index],
            this->property.get(index, source->wavelength_index),
            this->medium_property.get(index, source->wavelength_index),
            source
        );
    }

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const override {
        Sphere scatterer(
            this->diameter[index],
            this->property.get(index, source->wavelength_index),
            this->medium_property.get(index, source->wavelength_index),
            source
        );

        return std::make_unique<Sphere>(scatterer);
    }

    Sphere get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
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

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const override {
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
};


// Cylinder class inheriting from BaseSet
class CylinderSet : public ScattererSet {
public:
    std::vector<double> diameter;
    ScattererProperties property;
    MediumProperties medium_property;

    CylinderSet() = default;
    CylinderSet(const std::vector<double>& diameter, const ScattererProperties& property, const MediumProperties& medium_property, const bool is_sequential)
        : ScattererSet(is_sequential), diameter(diameter), property(property), medium_property(medium_property)
        {this->update_shape();}

    void update_shape() override {
        this->shape = {
            diameter.size(),
            property.size(),
            medium_property.size()
        };

        total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
    }

    void validate_sequential_data(const size_t expected_size) const override {
        // Check each vector's size and throw an error with the specific vector name if sizes don't match
        if (this->diameter.size() != expected_size)
            throw std::runtime_error("Error: Vector size mismatch in sequential computation. diameter has a different size than expected size.");

        if (this->property.size() != expected_size)
            throw std::runtime_error("Error: Vector size mismatch in sequential computation. property has a different size than expected size.");

        if (this->medium_property.size() != expected_size)
            throw std::runtime_error("Error: Vector size mismatch in sequential computation. medium_property has a different size than expected size.");
    }

    Cylinder get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
        return Cylinder(
            this->diameter[index],
            this->property.get(index, source->wavelength_index),
            this->medium_property.get(index, source->wavelength_index),
            source
        );
    }

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const override {
        Cylinder scatterer(
            this->diameter[index],
            this->property.get(index, source->wavelength_index),
            this->medium_property.get(index, source->wavelength_index),
            source
        );

        return std::make_unique<Cylinder>(scatterer);
    }

    Cylinder get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
        std::vector<size_t> indices = calculate_indices(flat_index);

        Cylinder scatterer(
            diameter[indices[0]],
            property.get(indices[1], source->wavelength_index),
            medium_property.get(indices[2], source->wavelength_index),
            source
        );

        scatterer.indices = indices;

        return scatterer;
    }


    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const override {
        std::vector<size_t> indices = calculate_indices(flat_index);

        Cylinder scatterer = Cylinder(
            diameter[indices[0]],
            property.get(indices[1], source->wavelength_index),
            medium_property.get(indices[2], source->wavelength_index),
            source
        );

        scatterer.indices = indices;

        return std::make_unique<Cylinder>(scatterer);
    }
};



// Core-shell class
class CoreShellSet : public ScattererSet {
public:
    std::vector<double> core_diameter;
    std::vector<double> shell_thickness;
    ScattererProperties core_property;
    ScattererProperties shell_property;
    MediumProperties medium_property;

    CoreShellSet() = default;

    CoreShellSet(
        const std::vector<double>& core_diameter,
        const std::vector<double>& shell_thickness,
        const ScattererProperties& core_property,
        const ScattererProperties& shell_property,
        const MediumProperties& medium_property,
        const bool is_sequential)
        : ScattererSet(is_sequential), core_diameter(core_diameter), shell_thickness(shell_thickness), core_property(core_property), shell_property(shell_property), medium_property(medium_property)
        {this->update_shape();}

    void update_shape() override {
        this->shape = {
            core_diameter.size(),
            shell_thickness.size(),
            core_property.size(),
            shell_property.size(),
            medium_property.size()
        };

        total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
    }

    CoreShell get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const {
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

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const override {
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

    void validate_sequential_data(const size_t expected_size) const override {
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

    CoreShell get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const {
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

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const override {
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
};
