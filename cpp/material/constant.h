#pragma once

#include "./base.h"

template<typename RefractiveIndexType>
class BaseConstant : public Base<RefractiveIndexType> {
public:
    RefractiveIndexType constant_refractive_index;

    BaseConstant() = default;

    explicit BaseConstant(const RefractiveIndexType& refractive_index)
        : constant_refractive_index(refractive_index) {
            this->refractive_index = refractive_index;
            this->is_initialized = true;
        }

    std::shared_ptr<Base<RefractiveIndexType>> clone() const override {
        return std::make_shared<BaseConstant<RefractiveIndexType>>(this->constant_refractive_index);
    }

protected:
    RefractiveIndexType compute_refractive_index(const double wavelength) const override {
        (void)wavelength;
        return constant_refractive_index;
    }
};


using ConstantMedium = BaseConstant<double>;
using ConstantMaterial = BaseConstant<complex128>;
