#pragma once

#include <cstdarg>    // For va_list, va_start, va_end
#include <cstdio>     // For printf, vprintf
#include "sphere/sphere.h"
#include "cylinder/cylinder.h"
#include "coreshell/coreshell.h"
#include "source/source.h"
#include "detector/detector.h"
#include "utils/numpy_interface.cpp"
#include "properties/properties.h"
#include "sets/detector.cpp"
#include "sets/scatterer.cpp"
#include "sets/source.cpp"

typedef std::complex<double> complex128;

#define DEFINE_SCATTERER_FUNCTION(scatterer, SCATTERER, dtype, name) \
    pybind11::array_t<dtype> get_##scatterer##_##name() const { \
        debug_printf("Computing %s_%s\n", #scatterer, #name); \
        return get_scatterer_data<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name); \
    } \
    pybind11::array_t<dtype> get_##scatterer##_##name##_sequential() const { \
        debug_printf("Computing %s_%s_sequential\n", #scatterer, #name); \
        return get_scatterer_data_sequential<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name); \
    }

#define DEFINE_SCATTERER_COEFFICIENT(scatterer, SCATTERER, name) \
    pybind11::array_t<double> get_##scatterer##_##name() const { \
        debug_printf("Computing %s_%s coefficient\n", #scatterer, #name); \
        return get_scatterer_data<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name##_abs); \
    } \
    pybind11::array_t<double> get_##scatterer##_##name##_sequential() const { \
        debug_printf("Computing %s_%s coefficient sequentially\n", #scatterer, #name); \
        return get_scatterer_data_sequential<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name##_abs); \
    }

#define DEFINE_SCATTERER_COUPLING(scatterer, SCATTERER) \
    pybind11::array_t<double> get_##scatterer##_coupling() const { \
        debug_printf("Computing %s_coupling\n", #scatterer); \
        return get_scatterer_coupling<SCATTERER::Set>(scatterer##Set); \
    } \
    pybind11::array_t<double> get_##scatterer##_coupling_sequential() const { \
        debug_printf("Computing %s_coupling_sequential\n", #scatterer); \
        return get_scatterer_coupling_sequential<SCATTERER::Set>(scatterer##Set); \
    }

class Experiment
{
    public:
        bool debug_mode = true;
        SPHERE::Set sphereSet;
        CYLINDER::Set cylinderSet;
        CORESHELL::Set coreshellSet;
        DETECTOR::Set detector_set;
        SOURCE::Set source_set;

        explicit Experiment(bool debug_mode = true) : debug_mode(debug_mode) {}

        // Setter methods
        void set_sphere(SPHERE::Set& set) { sphereSet = set; }
        void set_cylinder(CYLINDER::Set& set) { cylinderSet = set; }
        void set_coreshell(CORESHELL::Set& set) { coreshellSet = set; }
        void set_source(SOURCE::Set &set) { source_set = set; }
        void set_detector(DETECTOR::Set &set) { detector_set = set; }

        // Helper method for debugging without iostream.
        // It prints the given message only if debug_mode is true.
        void debug_printf(const char* format, ...) const {
            if (!debug_mode) return;
            va_list args;
            va_start(args, format);
            vprintf(format, args);
            va_end(args);
        }

        // Example: You can also add a class-wide debug_print method.
        void debug_print_state() const {
            if (!debug_mode) return;
            debug_printf("----- Experiment Debug Info -----\n");
            debug_printf("SourceSet total combinations: %zu\n", source_set.total_combinations);
            debug_printf("SphereSet total combinations: %zu\n", sphereSet.total_combinations);
            debug_printf("CylinderSet total combinations: %zu\n", cylinderSet.total_combinations);
            debug_printf("CoreshellSet total combinations: %zu\n", coreshellSet.total_combinations);
            debug_printf("DetectorSet total combinations: %zu\n", detector_set.total_combinations);
            debug_printf("SourceSet shape: ");
            for (size_t dim : source_set.shape) {
                debug_printf("%zu ", dim);
            }
            debug_printf("\n---------------------------------\n");
        }

        // Flattens a multi-index into a single index.
        template <typename... MultiIndexVectors>
        static size_t flatten_multi_index(const std::vector<size_t>& dimensions,const MultiIndexVectors&... multi_indices) {
            std::vector<size_t> multi_index = concatenate_vector(multi_indices...);
            size_t flatten_index = 0;
            size_t stride = 1;
            for (int i = static_cast<int>(dimensions.size()) - 1; i >= 0; --i) {
                flatten_index += multi_index[i] * stride;
                stride *= dimensions[i];
            }
            return flatten_index;
        }

        template<typename dtype, typename ScattererSet, typename Function>
        pybind11::array_t<dtype> get_scatterer_data(const ScattererSet& scatterer_set, Function function) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape;
            size_t total_iterations;

            if (detector_set.is_empty) {
                array_shape = concatenate_vector(source_set.shape, scatterer_set.shape);
                total_iterations = source_set.total_combinations * scatterer_set.total_combinations;
            }
            else {
                array_shape = concatenate_vector(source_set.shape, scatterer_set.shape, detector_set.shape);
                total_iterations = source_set.total_combinations * scatterer_set.total_combinations * detector_set.total_combinations;
            }

            debug_printf("get_scatterer_data: total_iterations = %zu\n", total_iterations);

            std::vector<dtype> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t flat_index = 0; flat_index < total_iterations; ++flat_index) {
                size_t idx; // Declare idx locally so each iteration has its own copy

                if (detector_set.is_empty) {
                    // 2D case: only source and scatterer
                    size_t i = flat_index / scatterer_set.total_combinations;
                    size_t j = flat_index % scatterer_set.total_combinations;
                    SOURCE::BaseSource source = source_set.get_source_by_index(i);
                    auto scatterer = scatterer_set.get_scatterer_by_index(j, source);
                    idx = flatten_multi_index(array_shape, source.indices, scatterer.indices);
                    dtype value = (scatterer.*function)();
                    output_array[idx] = value;
                } else {
                    // 3D case: source, scatterer, and detector
                    size_t i = flat_index / (scatterer_set.total_combinations * detector_set.total_combinations);
                    size_t j = (flat_index / detector_set.total_combinations) % scatterer_set.total_combinations;
                    size_t k = flat_index % detector_set.total_combinations;
                    SOURCE::BaseSource source = source_set.get_source_by_index(i);
                    auto scatterer = scatterer_set.get_scatterer_by_index(j, source);
                    DETECTOR::Detector detector = detector_set.get_detector_by_index(k);
                    idx = flatten_multi_index(array_shape, source.indices, scatterer.indices, detector.indices);
                    dtype value = (scatterer.*function)();
                    output_array[idx] = value;
                }
            }

            debug_printf("get_scatterer_data: finished computation\n");
            return _vector_to_numpy(output_array, array_shape);
        }

        template<typename dtype, typename ScattererSet, typename Function>
        pybind11::array_t<dtype> get_scatterer_data_sequential(const ScattererSet& scatterer_set, Function function) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape = {source_set.wavelength.size()};
            size_t full_size = source_set.wavelength.size();
            scatterer_set.validate_sequential_data(full_size);
            source_set.validate_sequential_data(full_size);
            debug_printf("get_scatterer_data_sequential: full_size = %zu\n", full_size);

            std::vector<dtype> output_array(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {
                SOURCE::BaseSource source = source_set.get_source_by_index_sequential(idx);
                auto scatterer = scatterer_set.get_scatterer_by_index_sequential(idx, source);
                output_array[idx] = (scatterer.*function)();
            }
            debug_printf("get_scatterer_data_sequential: finished computation\n");
            return _vector_to_numpy(output_array, {full_size});
        }

        template<typename ScattererSet>
        pybind11::array_t<double> get_scatterer_coupling(const ScattererSet& scatterer_set) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape = concatenate_vector(source_set.shape, scatterer_set.shape, detector_set.shape);
            size_t total_iterations = source_set.total_combinations * scatterer_set.total_combinations * detector_set.total_combinations;
            debug_printf("get_scatterer_coupling: total_iterations = %zu\n", total_iterations);

            std::vector<double> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t idx_flat = 0; idx_flat < total_iterations; ++idx_flat) {
                size_t i = idx_flat / (scatterer_set.total_combinations * detector_set.total_combinations);
                size_t j = (idx_flat / detector_set.total_combinations) % scatterer_set.total_combinations;
                size_t k = idx_flat % detector_set.total_combinations;
                SOURCE::BaseSource source = source_set.get_source_by_index(i);
                auto scatterer = scatterer_set.get_scatterer_by_index(j, source);
                DETECTOR::Detector detector = detector_set.get_detector_by_index(k);
                detector.medium_refractive_index = scatterer.medium_refractive_index;
                size_t idx = flatten_multi_index(array_shape, source.indices, scatterer.indices, detector.indices);
                output_array[idx] = detector.get_coupling(scatterer);
            }
            debug_printf("get_scatterer_coupling: finished computation\n");
            return _vector_to_numpy(output_array, array_shape);
        }

        template<typename ScattererSet>
        pybind11::array_t<double> get_scatterer_coupling_sequential(const ScattererSet& scatterer_set) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape = {source_set.wavelength.size()};
            size_t full_size = source_set.wavelength.size();
            scatterer_set.validate_sequential_data(full_size);
            source_set.validate_sequential_data(full_size);
            detector_set.validate_sequential_data(full_size);
            debug_printf("get_scatterer_coupling_sequential: full_size = %zu\n", full_size);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {
                SOURCE::BaseSource source = source_set.get_source_by_index_sequential(idx);
                auto scatterer = scatterer_set.get_scatterer_by_index_sequential(idx, source);
                DETECTOR::Detector detector = detector_set.get_detector_by_index_sequential(idx);
                output_array[idx] = detector.get_coupling(scatterer);
            }
            debug_printf("get_scatterer_coupling_sequential: finished computation\n");
            return _vector_to_numpy(output_array, {full_size});
        }

        //--------------------------------------SPHERE------------------------------------
        DEFINE_SCATTERER_COEFFICIENT(sphere, SPHERE, a1)
        DEFINE_SCATTERER_COEFFICIENT(sphere, SPHERE, a2)
        DEFINE_SCATTERER_COEFFICIENT(sphere, SPHERE, a3)
        DEFINE_SCATTERER_COEFFICIENT(sphere, SPHERE, b1)
        DEFINE_SCATTERER_COEFFICIENT(sphere, SPHERE, b2)
        DEFINE_SCATTERER_COEFFICIENT(sphere, SPHERE, b3)

        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Qsca)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Qext)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Qabs)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Qpr)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Qback)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Qforward)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Qratio)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Csca)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Cext)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Cabs)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Cpr)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Cback)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Cratio)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, Cforward)
        DEFINE_SCATTERER_FUNCTION(sphere, SPHERE, double, g)

        DEFINE_SCATTERER_COUPLING(sphere, SPHERE)

        //--------------------------------------CYLINDER------------------------------------
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, a11)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, a12)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, a13)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, a21)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, a22)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, a23)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, b11)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, b12)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, b13)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, b21)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, b22)
        DEFINE_SCATTERER_COEFFICIENT(cylinder, CYLINDER, b23)

        DEFINE_SCATTERER_FUNCTION(cylinder, CYLINDER, double, Qsca)
        DEFINE_SCATTERER_FUNCTION(cylinder, CYLINDER, double, Qext)
        DEFINE_SCATTERER_FUNCTION(cylinder, CYLINDER, double, Qabs)
        DEFINE_SCATTERER_FUNCTION(cylinder, CYLINDER, double, Csca)
        DEFINE_SCATTERER_FUNCTION(cylinder, CYLINDER, double, Cext)
        DEFINE_SCATTERER_FUNCTION(cylinder, CYLINDER, double, Cabs)
        DEFINE_SCATTERER_FUNCTION(cylinder, CYLINDER, double, g)

        DEFINE_SCATTERER_COUPLING(cylinder, CYLINDER)

        //--------------------------------------CORESHELL------------------------------------
        DEFINE_SCATTERER_COEFFICIENT(coreshell, CORESHELL, a1)
        DEFINE_SCATTERER_COEFFICIENT(coreshell, CORESHELL, a2)
        DEFINE_SCATTERER_COEFFICIENT(coreshell, CORESHELL, a3)
        DEFINE_SCATTERER_COEFFICIENT(coreshell, CORESHELL, b1)
        DEFINE_SCATTERER_COEFFICIENT(coreshell, CORESHELL, b2)
        DEFINE_SCATTERER_COEFFICIENT(coreshell, CORESHELL, b3)

        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Qsca)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Qext)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Qabs)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Qpr)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Qback)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Qforward)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Qratio)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Csca)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Cext)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Cabs)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Cpr)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Cback)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Cratio)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, Cforward)
        DEFINE_SCATTERER_FUNCTION(coreshell, CORESHELL, double, g)

        DEFINE_SCATTERER_COUPLING(coreshell, CORESHELL)
};
