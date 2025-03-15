#pragma once

#include "sphere/sphere.h"
#include "cylinder/cylinder.h"
#include "coreshell/coreshell.h"
#include "detector/detector.h"
#include "source/source.h"


// Base class to reduce redundancy
class BaseSet{
    public:
        bool is_sequential;
        std::vector<size_t> shape = {1};
        size_t current_index = 0;
        size_t total_combinations = 1;
        bool is_empty = true;
        virtual void update_shape() {};

        BaseSet() = default;
        BaseSet(bool is_sequential);
        virtual ~BaseSet() = default;

        // Calculate the multi-dimensional indices for the current index
        std::vector<size_t> calculate_indices(size_t flat_index) const;

        template <typename Container> void check_size(const Container& container, size_t expected_size, const std::string& name) const;
};
