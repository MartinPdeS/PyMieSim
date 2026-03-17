#pragma once

#include <single/scatterer/coreshell/coreshell.h>
#include <experiment/scatterer_set/base_set.h>
#include <experiment/material_set/material_set.h>


class CoreShellSet : public ScattererSet {
public:
    inline static const std::vector<std::string> attributes = {
        "core_diameter",
        "shell_thickness",
        "core_material",
        "shell_material",
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

    std::vector<double> core_diameter;
    std::vector<double> shell_thickness;
    MaterialSet core_material;
    MaterialSet shell_material;
    MediumSet medium;

    CoreShellSet() = default;

    CoreShellSet(
        const std::vector<double>& core_diameter,
        const std::vector<double>& shell_thickness,
        const MaterialSet& core_material,
        const MaterialSet& shell_material,
        const MediumSet& medium,
        const bool is_sequential)
        : ScattererSet(is_sequential), core_diameter(core_diameter), shell_thickness(shell_thickness), core_material(core_material), shell_material(shell_material), medium(medium)
        {this->update_shape();}

    void update_shape() override;

    std::shared_ptr<BaseScatterer> get_scatterer_by_index(const size_t flat_index) const override;

    std::shared_ptr<BaseScatterer> get_scatterer_by_index_sequential(const size_t index) const override;
};
