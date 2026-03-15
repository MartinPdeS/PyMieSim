#pragma once

#include <memory>
#include "experiment/base_set.h"
#include "single/scatterer/base_scatterer/base_scatterer.h"
#include <material/material.h>


class ScattererSet: public BaseSet
{
public:
    ScattererSet() = default;
    ScattererSet(const bool is_sequential) : BaseSet(is_sequential) {}

    virtual void validate_sequential_data(const size_t expected_size) const = 0;
    virtual std::shared_ptr<BaseScatterer> get_scatterer_by_index_sequential(const size_t index) const = 0;
    virtual std::shared_ptr<BaseScatterer> get_scatterer_by_index(const size_t flat_index) const = 0;

};

