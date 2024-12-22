#pragma once

#include "single/includes/sphere.cpp"
#include "single/includes/cylinder.cpp"
#include "single/includes/coreshell.cpp"
#include "single/includes/sources.cpp"
#include "single/includes/detectors.cpp"

#include "utils/numpy_interface.cpp"
#include "experiment/includes/scatterer_properties.cpp"
#include "experiment/includes/sets.cpp"

typedef std::complex<double> complex128;


#define DEFINE_SCATTERER_FUNCTION(scatterer, SCATTERER, dtype, name) \
    pybind11::array_t<dtype> get_##scatterer##_##name() const { return get_scatterer_data<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name); } \
    pybind11::array_t<dtype> get_##scatterer##_##name##_sequential() const { return get_scatterer_data_sequential<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name); }

#define DEFINE_SCATTERER_COEFFICIENT(scatterer, SCATTERER, name) \
    pybind11::array_t<double> get_##scatterer##_##name() const {return get_scatterer_data<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name##_abs);} \
    pybind11::array_t<double> get_##scatterer##_##name##_sequential() const {return get_scatterer_data_sequential<double, SCATTERER::Set>(scatterer##Set, &SCATTERER::Scatterer::get_##name##_abs);}

#define DEFINE_SCATTERER_COUPLING(scatterer, SCATTERER) \
    pybind11::array_t<double> get_##scatterer##_coupling() const { return get_scatterer_coupling<SCATTERER::Set>(scatterer##Set); } \
    pybind11::array_t<double> get_##scatterer##_coupling_sequential() const { return get_scatterer_coupling_sequential<SCATTERER::Set>(scatterer##Set); }

class Experiment
{
    public:
        SPHERE::Set sphereSet;
        CYLINDER::Set cylinderSet;
        CORESHELL::Set coreshellSet;
        DETECTOR::Set detectorSet;
        SOURCE::Set sourceSet;

        Experiment() = default;

        void set_sphere(SPHERE::Set& set) { sphereSet = set; }
        void set_cylinder(CYLINDER::Set& set) { cylinderSet = set; }
        void set_coreshell(CORESHELL::Set& set) { coreshellSet = set; }
        void set_source(SOURCE::Set &set) { sourceSet = set; }
        void set_detector(DETECTOR::Set &set) { detectorSet = set; }

        template <typename... MultiIndexVectors>
        static size_t flatten_multi_index(const std::vector<size_t>& dimensions, const MultiIndexVectors&... multi_indices) {

            // Concatenate all the multi_index vectors into one
            std::vector<size_t> multi_index = concatenate_vector(multi_indices...);
            size_t flatten_index = 0;
            size_t stride = 1;

            const size_t* multi_index_ptr = multi_index.data();
            const size_t* dimensions_ptr = dimensions.data();

            // Iterate from the last dimension to the first
            for (int i = dimensions.size() - 1; i >= 0; --i) {
                flatten_index += multi_index_ptr[i] * stride;
                stride *= dimensions_ptr[i];
            }

            return flatten_index;
        }

        template<typename dtype, typename ScattererSet, typename Function>
        pybind11::array_t<dtype> get_scatterer_data(const ScattererSet& scattererSet, Function function) const {
            std::vector<size_t> array_shape = concatenate_vector(sourceSet.shape, scattererSet.shape);

            size_t total_iterations = sourceSet.total_combinations * scattererSet.total_combinations;

            std::vector<dtype> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t flat_index = 0; flat_index < total_iterations; ++flat_index) {
                // Map flat_index to i and j
                size_t i = flat_index / scattererSet.total_combinations;
                size_t j = flat_index % scattererSet.total_combinations;

                // Retrieve source and scatterer for the current combination of i and j
                SOURCE::BaseSource source = sourceSet.get_source_by_index(i);

                auto scatterer = scattererSet.get_scatterer_by_index(j, source);

                // Flatten multi-dimensional index for output
                size_t idx = flatten_multi_index(array_shape, source.indices, scatterer.indices);

                // Perform the computation and store the result
                dtype value = (scatterer.*function)();
                output_array[idx] = value;
            }

            return _vector_to_numpy(output_array, array_shape);
        }

        template<typename dtype, typename ScattererSet, typename Function>
        pybind11::array_t<dtype> get_scatterer_data_sequential(const ScattererSet& scattererSet, Function function) const {

            std::vector<size_t> array_shape = {sourceSet.wavelength.size()};

            size_t full_size = sourceSet.wavelength.size();

            scattererSet.validate_sequential_data(full_size);
            sourceSet.validate_sequential_data(full_size);

            std::vector<dtype> output_array(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {
                SOURCE::BaseSource source = sourceSet.get_source_by_index_sequential(idx);

                auto scatterer = scattererSet.get_scatterer_by_index_sequential(idx, source);

                output_array[idx] = (scatterer.*function)();
            }

            return _vector_to_numpy(output_array, {full_size});
        }


        template<typename ScattererSet>
        pybind11::array_t<double> get_scatterer_coupling(const ScattererSet& scattererSet) const {
            std::vector<size_t> array_shape = concatenate_vector(sourceSet.shape, scattererSet.shape, detectorSet.shape);

            // Calculate total iterations for the single loop
            size_t total_iterations = sourceSet.total_combinations * scattererSet.total_combinations * detectorSet.total_combinations;

            std::vector<double> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t idx_flat = 0; idx_flat < total_iterations; ++idx_flat) {
                // Calculate i, j, k from idx_flat
                size_t i = idx_flat / (scattererSet.total_combinations * detectorSet.total_combinations);
                size_t j = (idx_flat / detectorSet.total_combinations) % scattererSet.total_combinations;
                size_t k = idx_flat % detectorSet.total_combinations;

                // Retrieve the source, scatterer, and detector based on the calculated indices
                SOURCE::BaseSource source = sourceSet.get_source_by_index(i);
                auto scatterer = scattererSet.get_scatterer_by_index(j, source);
                DETECTOR::Detector detector = detectorSet.get_detector_by_index(k);

                // Flatten the multi-index into a single index for output_array
                size_t idx = flatten_multi_index(array_shape, source.indices, scatterer.indices, detector.indices);

                // Perform the operation
                output_array[idx] = detector.get_coupling(scatterer);
            }

            return _vector_to_numpy(output_array, array_shape);
        }

        template<typename ScattererSet>
        pybind11::array_t<double> get_scatterer_coupling_sequential(const ScattererSet& scattererSet) const {
            std::vector<size_t> array_shape = {sourceSet.wavelength.size()};

            size_t full_size = sourceSet.wavelength.size();

            std::vector<double> output_array(full_size);

            scattererSet.validate_sequential_data(full_size);
            sourceSet.validate_sequential_data(full_size);
            detectorSet.validate_sequential_data(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {
                SOURCE::BaseSource source = sourceSet.get_source_by_index_sequential(idx);

                auto scatterer = scattererSet.get_scatterer_by_index_sequential(idx, source);

                DETECTOR::Detector detector = detectorSet.get_detector_by_index_sequential(idx);

                output_array[idx] = detector.get_coupling(scatterer);
            }
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


