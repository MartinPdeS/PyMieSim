#pragma once


#include <experiment/sets/scatterer_set/base.h>



// InfiniteCylinder class inheriting from BaseSet
class InfiniteCylinderSet : public ScattererSet {
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
    InfiniteCylinderSet(
        const std::vector<double>& diameter,
        const ScattererProperties& property,
        const MediumProperties& medium_property,
        const bool is_sequential
    )
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
