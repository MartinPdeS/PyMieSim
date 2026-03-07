#pragma once

#include <memory>
#include "experiment/sets/base_set.h"
#include "single/scatterer/base_scatterer/base_scatterer.h"
#include "experiment/sets/properties_set/properties_set.h"
#include <dispersive_material/material.h>


class ScattererSet: public BaseSet
{
public:
    ScattererSet() = default;
    ScattererSet(const bool is_sequential) : BaseSet(is_sequential) {}

    virtual void validate_sequential_data(const size_t expected_size) const = 0;
    virtual std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index_sequential(const size_t index, std::shared_ptr<BaseSource> source) const = 0;
    virtual std::unique_ptr<BaseScatterer> get_scatterer_ptr_by_index(const size_t flat_index, std::shared_ptr<BaseSource> source) const = 0;

};

