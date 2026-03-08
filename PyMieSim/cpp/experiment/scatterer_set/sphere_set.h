#pragma once

#include <experiment/scatterer_set/base_set.h>
#include <experiment/material_set/material_set.h>
#include <experiment/material_set/medium_set.h>


// Sphere class inheriting from BaseSet
class SphereSet : public ScattererSet
{
public:
    inline static const std::vector<std::string> attributes = {
        "diameter",
        "property",
        "medium_property"
    };
    inline static const std::vector<std::string> available_measure_list = {
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
    MaterialSet property;
    MediumSet medium_property;

    SphereSet() = default;

    SphereSet(
        const std::vector<double>& diameter,
        const MaterialSet& property,
        const MediumSet& medium_property,
        const bool is_sequential
    )
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