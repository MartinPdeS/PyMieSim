#pragma once

#include <memory>
#include "experiment/sets/base_set.h"
#include "single/scatterer/base_scatterer/base_scatterer.h"
#include "experiment/sets/properties_set/properties_set.h"

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
    std::vector<std::string> attributes = {
        "diameter",
        "property",
        "medium_property"
    };
    std::vector<std::string> available_measure_list = {
        "Qsca",
        "Qext",
        "Qabs",
        "Qratio",
        "Qforward",
        "Qback",
        "Qpr",
        "Csca",
        "Cext",
        "Cabs",
        "Cratio",
        "Cforward",
        "Cback",
        "Cpr",
        "a1",
        "a2",
        "a3",
        "b1",
        "b2",
        "b3",
        "g",
        "coupling",
    };


    std::vector<double> diameter;
    ScattererProperties property;
    MediumProperties medium_property;

    SphereSet() = default;

    SphereSet(const std::vector<double>& diameter, const ScattererProperties& property, const MediumProperties& medium_property, const bool is_sequential)
        : ScattererSet(is_sequential), diameter(diameter), property(property), medium_property(medium_property)
        {
            this->update_shape();
        }

    void update_shape() override;

    void validate_sequential_data(const size_t expected_size) const override;

    Sphere get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const;

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const override;

    Sphere get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const;

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const override;
};


// InfiniteCylinder class inheriting from BaseSet
class InfiniteCylinderSet : public ScattererSet {
public:
    std::vector<std::string> attributes = {
        "diameter",
        "property",
        "medium_property"
    };

    std::vector<std::string> available_measure_list = {
        "Qsca",
        "Qext",
        "Qabs",
        "Csca",
        "Cext",
        "Cabs",
        "a11",
        "a21",
        "a12",
        "a22",
        "a13",
        "a23",
        "b11",
        "b21",
        "b12",
        "b22",
        "b13",
        "b23",
        "coupling",
    };
    std::vector<double> diameter;
    ScattererProperties property;
    MediumProperties medium_property;

    InfiniteCylinderSet() = default;
    InfiniteCylinderSet(const std::vector<double>& diameter, const ScattererProperties& property, const MediumProperties& medium_property, const bool is_sequential)
        : ScattererSet(is_sequential), diameter(diameter), property(property), medium_property(medium_property)
        {
            this->update_shape();
        }

    void update_shape() override;

    void validate_sequential_data(const size_t expected_size) const override;

    InfiniteCylinder get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const;

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const override;

    InfiniteCylinder get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const;


    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const override;
};



// Core-shell class
class CoreShellSet : public ScattererSet {
public:
    std::vector<std::string> attributes = {
        "core_diameter",
        "shell_thickness",
        "core_property",
        "shell_property",
        "medium_property"
    };

    std::vector<std::string> available_measure_list = {
        "Qsca",
        "Qext",
        "Qabs",
        "Qratio",
        "Qforward",
        "Qback",
        "Qpr",
        "Csca",
        "Cext",
        "Cabs",
        "Cratio",
        "Cforward",
        "Cback",
        "Cpr",
        "a1",
        "a2",
        "a3",
        "b1",
        "b2",
        "b3",
        "g",
        "coupling",
    };

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

    void update_shape() override;

    CoreShell get_scatterer_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const;

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const override;

    void validate_sequential_data(const size_t expected_size) const override;

    CoreShell get_scatterer_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const;

    std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const override;
};
