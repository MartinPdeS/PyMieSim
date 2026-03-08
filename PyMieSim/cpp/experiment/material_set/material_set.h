#pragma once

#include <single/material/material.h>

class MaterialSet
{
public:
    std::vector<std::shared_ptr<BaseMaterial>> materials;
    MaterialSet() = default;

    MaterialSet(
        const std::vector<std::shared_ptr<BaseMaterial>>& materials
    ): materials(materials) {}

    MaterialSet(
        const std::vector<complex128>& refractive_indices
    )
    {
        for (const auto& ri : refractive_indices) {
            auto material = std::make_shared<ConstantMaterial>(ri);
            this->materials.push_back(material);
        }
    }

    size_t size() const {
        return this->materials.size();
    }

    std::shared_ptr<BaseMaterial> operator[](size_t index) const {
        return this->materials[index];
    }
};

