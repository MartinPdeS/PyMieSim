#pragma once

#include <vector>
#include <memory>
#include <stdexcept>

#include <material/material.h>

class MaterialSet
{
public:
    std::vector<std::shared_ptr<BaseMaterial>> materials;

    MaterialSet() = default;

    MaterialSet(
        const std::vector<std::shared_ptr<BaseMaterial>>& materials
    )
        : materials(materials)
    {}

    MaterialSet(
        const std::vector<complex128>& refractive_indices
    )
    {
        this->materials.reserve(refractive_indices.size());

        for (const auto& refractive_index : refractive_indices) {
            this->materials.push_back(
                std::make_shared<ConstantMaterial>(refractive_index)
            );
        }
    }

    size_t size() const {
        return this->materials.size();
    }

    bool empty() const {
        return this->materials.empty();
    }

    std::shared_ptr<BaseMaterial> operator[](size_t index) const {
        return this->materials[index];
    }
};

class MediumSet
{
public:
    std::vector<std::shared_ptr<BaseMedium>> mediums;

    MediumSet() = default;

    MediumSet(
        const std::vector<std::shared_ptr<BaseMedium>>& mediums
    )
        : mediums(mediums)
    {}

    MediumSet(
        const std::vector<double>& refractive_indices
    )
    {
        this->mediums.reserve(refractive_indices.size());

        for (const auto& refractive_index : refractive_indices) {
            this->mediums.push_back(
                std::make_shared<ConstantMedium>(refractive_index)
            );
        }
    }

    size_t size() const {
        return this->mediums.size();
    }

    bool empty() const {
        return this->mediums.empty();
    }

    std::shared_ptr<BaseMedium> operator[](size_t index) const {
        return this->mediums[index];
    }
};