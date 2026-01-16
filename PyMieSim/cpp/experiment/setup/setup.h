#pragma once

#include <functional>
#include <memory>
#include <cstdarg>    // For va_list, va_start, va_end
#include <cstdio>     // For printf, vprintf
#include "single/scatterer/sphere/sphere.h"
#include "single/scatterer/cylinder/cylinder.h"
#include "single/scatterer/coreshell/coreshell.h"
#include "single/source/source.h"
#include "single/detector/detector.h"
#include "experiment/sets/sets.cpp"

class Setup
{
    public:
        bool debug_mode = false;

        explicit Setup(bool debug_mode = false) : debug_mode(debug_mode) {}


        // Helper method for debugging without iostream.
        // It prints the given message only if debug_mode is true.
        void debug_printf(const char* format, ...) const {
            if (!debug_mode) return;
            va_list args;
            va_start(args, format);
            vprintf(format, args);
            va_end(args);
        }

        /**
         * @brief Prints the current state of the experiment for debugging purposes.
         * @param scatterer_set The set of scatterers involved in the experiment.
         * @param source_set The set of sources involved in the experiment.
         * @param detector_set The set of detectors involved in the experiment.
        */
        void debug_print_state(const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const DetectorSet& detector_set) const;

        /**
         * @brief Flattens a multi-dimensional index into a single index (one input vector).
         */
        std::size_t flatten_multi_index(
            const std::vector<std::size_t>& dimensions,
            const std::vector<std::size_t>& multi_index
        ) const
        {
            if (multi_index.size() != dimensions.size()) {
                throw std::invalid_argument("flatten_multi_index: multi_index.size() must match dimensions.size().");
            }

            std::size_t flattened_index = 0;
            std::size_t stride = 1;

            for (std::size_t i = dimensions.size(); i-- > 0;) {
                flattened_index += multi_index[i] * stride;
                stride *= dimensions[i];
            }

            return flattened_index;
        }

        /**
         * @brief Flattens a multi-dimensional index into a single index (two input vectors).
         */
        std::size_t flatten_multi_index(
            const std::vector<std::size_t>& dimensions,
            const std::vector<std::size_t>& first_multi_index,
            const std::vector<std::size_t>& second_multi_index
        ) const
        {
            const std::vector<std::size_t> multi_index =
                this->concatenate_vector(first_multi_index, second_multi_index);

            return this->flatten_multi_index(dimensions, multi_index);
        }

        /**
         * @brief Flattens a multi-dimensional index into a single index (three input vectors).
         */
        std::size_t flatten_multi_index(
            const std::vector<std::size_t>& dimensions,
            const std::vector<std::size_t>& first_multi_index,
            const std::vector<std::size_t>& second_multi_index,
            const std::vector<std::size_t>& third_multi_index
        ) const
        {
            const std::vector<std::size_t> multi_index =
                this->concatenate_vector(first_multi_index, second_multi_index, third_multi_index);
            return this->flatten_multi_index(dimensions, multi_index);
        }

        /**
         * @brief Concatenates two vectors into a single vector.
         */
        std::vector<std::size_t> concatenate_vector(const std::vector<std::size_t> first_vector, const std::vector<std::size_t>& second_vector) const
        {
            std::vector<std::size_t> output_vector;
            output_vector.reserve(first_vector.size() + second_vector.size());

            output_vector.insert(output_vector.end(), first_vector.begin(), first_vector.end());
            output_vector.insert(output_vector.end(), second_vector.begin(), second_vector.end());

            return output_vector;
        }

        /**
         * @brief Concatenates three vectors into a single vector.
         */
        std::vector<std::size_t> concatenate_vector(const std::vector<std::size_t>& first_vector, const std::vector<std::size_t>& second_vector, const std::vector<std::size_t>& third_vector) const
        {
            std::vector<std::size_t> output_vector;
            output_vector.reserve(first_vector.size() + second_vector.size() + third_vector.size());

            output_vector.insert(output_vector.end(), first_vector.begin(), first_vector.end());
            output_vector.insert(output_vector.end(), second_vector.begin(), second_vector.end());
            output_vector.insert(output_vector.end(), third_vector.begin(), third_vector.end());

            return output_vector;
        }


        /**
         * @brief Generic method to get scatterer data by invoking a member function.
         * @tparam function Pointer to the member function of BaseScatterer to invoke.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A tuple containing a numpy array of the requested data and the shape of the array.
         */
        template<double (BaseScatterer::*function)() const>
        std::tuple<std::vector<double>, std::vector<std::size_t>>
        get_data(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const {

            if (debug_mode)
                this->debug_print_state(scatterer_set, source_set, detector_set);

            std::vector<std::size_t> array_shape;
            size_t total_iterations;

            if (detector_set.is_empty) {
                array_shape = this->concatenate_vector(source_set.shape, scatterer_set.shape);
                total_iterations = source_set.total_combinations * scatterer_set.total_combinations;
            }
            else {
                array_shape = this->concatenate_vector(source_set.shape, scatterer_set.shape, detector_set.shape);
                total_iterations = source_set.total_combinations * scatterer_set.total_combinations * detector_set.total_combinations;
            }

            debug_printf("get_scatterer_data: total_iterations = %zu\n", total_iterations);

            std::vector<double> output_array(total_iterations);

            #pragma omp parallel for
            for (long long flat_index = 0; flat_index < static_cast<long long>(total_iterations); ++flat_index) {
                size_t idx; // Declare idx locally so each iteration has its own copy

                if (detector_set.is_empty) {
                    // 2D case: only source and scatterer
                    size_t i = flat_index / scatterer_set.total_combinations;
                    size_t j = flat_index % scatterer_set.total_combinations;
                    std::shared_ptr<BaseSource> source_ptr = source_set.get_source_by_index(i);

                    std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index(j, source_ptr);

                    const std::vector<std::size_t> multi_index =
                        this->concatenate_vector(source_ptr->indices, source_ptr->indices);

                    idx = this->flatten_multi_index(array_shape, source_ptr->indices, scatterer_ptr->indices);
                    output_array[idx] = std::invoke(function, *scatterer_ptr);
                } else {
                    // 3D case: source, scatterer, and detector
                    long long i = flat_index / (scatterer_set.total_combinations * detector_set.total_combinations);
                    long long j = (flat_index / detector_set.total_combinations) % scatterer_set.total_combinations;
                    long long k = flat_index % detector_set.total_combinations;
                    std::shared_ptr<BaseSource> source_ptr = source_set.get_source_by_index(i);

                    std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index(j, source_ptr);

                    Detector detector = detector_set.get_detector_by_index(k);
                    idx = this->flatten_multi_index(array_shape, source_ptr->indices, scatterer_ptr->indices, detector.indices);

                    output_array[idx] = std::invoke(function, *scatterer_ptr);
                }
            }

            debug_printf("get_scatterer_data: finished computation\n");
            return std::make_tuple(std::move(output_array), std::move(array_shape));
        }


        /**
         * @brief Generic method to get scatterer data sequentially by invoking a member function.
         * @tparam function Pointer to the member function of BaseScatterer to invoke.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A numpy array containing the requested data.
         */
        template<double (BaseScatterer::*function)() const >
        std::vector<double>
        get_data_sequential(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const {

            if (debug_mode)
                this->debug_print_state(scatterer_set, source_set, detector_set);

            std::vector<std::size_t> array_shape = {source_set.wavelength.size()};
            size_t full_size = source_set.wavelength.size();
            scatterer_set.validate_sequential_data(full_size);
            source_set.validate_sequential_data(full_size);
            debug_printf("get_scatterer_data_sequential: full_size = %zu\n", full_size);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for
            for (long long idx = 0; idx < static_cast<long long>(full_size); ++idx) {
                std::shared_ptr<BaseSource> source = source_set.get_source_by_index_sequential(idx);

                std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index_sequential(idx, source);

                output_array[idx] = std::invoke(function, *scatterer_ptr);

            }
            debug_printf("get_scatterer_data_sequential: finished computation\n");
            return output_array;
        }

        /**
         * @brief Computes the coupling coefficient for given scatterers, sources, and detectors.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A tuple containing a numpy array of coupling coefficients and the shape of the array.
         */
        std::tuple<std::vector<double>, std::vector<std::size_t>>
        get_coupling(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const;

        /**
         * @brief Computes the coupling coefficient sequentially for given scatterers, sources, and detectors.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A numpy array of coupling coefficients.
         */
        std::vector<double>
        get_coupling_sequential(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const;

        /**
         * @brief Computes the far-field patterns for given scatterers, sources, and a Fibonacci mesh.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param mesh The Fibonacci mesh for far-field sampling.
         * @param distance The distance at which to compute the far-fields (default is 1).
         * @return A tuple containing a numpy array of far-field patterns and the shape of the array.
         */
        std::tuple<std::vector<complex128>, std::vector<std::size_t>>
        get_farfields(const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const FibonacciMesh& mesh, const double distance = 1) const;

};
