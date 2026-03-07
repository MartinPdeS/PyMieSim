#pragma once


#include <experiment/sets/scatterer_set/base.h>

// Core-shell class
class CoreShellSet : public ScattererSet {
public:
    inline static const std::vector<std::string> attributes = {
        "core_diameter",
        "shell_thickness",
        "core_property",
        "shell_property",
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
