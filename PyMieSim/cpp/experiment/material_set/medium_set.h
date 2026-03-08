#pragma once

#include <single/material/medium.h>

class MediumSet
{
public:
    std::vector<std::shared_ptr<BaseMedium>> mediums;
    MediumSet() = default;

    MediumSet(
        const std::vector<std::shared_ptr<BaseMedium>>& mediums
    ): mediums(mediums) {}

    MediumSet(
        const std::vector<double>& refractive_indices
    )
    {
        for (const auto& ri : refractive_indices) {
            auto medium = std::make_shared<ConstantMedium>(ri);
            this->mediums.push_back(medium);
        }
    }

    size_t size() const {
        return this->mediums.size();
    }

    std::shared_ptr<BaseMedium> operator[](size_t index) const {
        return this->mediums[index];
    }
};

