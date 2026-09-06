#pragma once

#include <single/scatterer/sphere/sphere.h>
#include <experiment/scatterer_set/base_set.h>
#include <experiment/material_set/material_set.h>

class SphereSet : public ScattererSet
{
public:
    inline static const std::vector<std::string> attributes = {
        "diameter",
        "material",
        "medium"
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
    MaterialSet material;
    MediumSet medium;

    SphereSet() = default;

    SphereSet(
        const std::vector<double>& diameter,
        const MaterialSet& material,
        const MediumSet& medium,
        const bool is_sequential
    )
        : ScattererSet(is_sequential), diameter(diameter), material(material), medium(medium)
    {
        this->update_shape();
    }

    void update_shape() override;

    std::shared_ptr<BaseScatterer> get_scatterer_by_index_sequential(const size_t index) const override;

    std::shared_ptr<BaseScatterer> get_scatterer_by_index(const size_t flat_index) const override;

};