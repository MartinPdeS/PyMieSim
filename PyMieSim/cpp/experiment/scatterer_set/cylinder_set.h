#pragma once


#include <experiment/scatterer_set/base_set.h>
#include <experiment/material_set/material_set.h>


// InfiniteCylinder class inheriting from BaseSet
class InfiniteCylinderSet : public ScattererSet {
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
    MaterialSet material;
    MediumSet medium;

    InfiniteCylinderSet() = default;
    InfiniteCylinderSet(
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

    void validate_sequential_data(const size_t expected_size) const override;

    std::shared_ptr<BaseScatterer> get_scatterer_by_index_sequential(const size_t index) const;

    std::shared_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index) const override;

    std::shared_ptr<BaseScatterer> get_scatterer_by_index(const size_t flat_index) const;

    std::shared_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index) const override;
};
