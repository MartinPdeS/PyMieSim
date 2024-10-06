#pragma once

#include "single/includes/sphere.cpp"

template<typename dtype, typename Function>
pybind11::array_t<dtype> Experiment::get_sphere_data(Function function) const
{
    using namespace SPHERE;
    std::vector<size_t> array_shape = concatenate_vector(sourceSet.shape, sphereSet.shape);
    size_t full_size = get_vector_sigma(array_shape);
    std::vector<dtype> output_array(full_size);

    #pragma omp parallel for
    for (size_t i = 0; i < sourceSet.total_combinations; ++i)
    {
        SOURCE::BaseSource source = sourceSet.get_source_by_index(i);

        #pragma omp parallel for
        for (size_t j = 0; j < sphereSet.total_combinations; ++j)
        {
            Scatterer scatterer = sphereSet.get_scatterer_by_index(j, source);
            size_t idx = flatten_multi_index(source.indices, scatterer.indices, array_shape);

            output_array[idx] = (scatterer.*function)();
        }
    }

    return vector_to_numpy(output_array, array_shape);
}


pybind11::array_t<double> Experiment::get_sphere_coupling() const
{
    using namespace SPHERE;
    using namespace SOURCE;
    using namespace DETECTOR;

    std::vector<size_t> array_shape = concatenate_vector(sourceSet.shape, sphereSet.shape, detectorSet.shape);
    size_t full_size = get_vector_sigma(array_shape);
    std::vector<double> output_array(full_size);

    #pragma omp parallel for
    for (size_t i = 0; i < sourceSet.total_combinations; ++i)
    {
        BaseSource source = sourceSet.get_source_by_index(i);

        #pragma omp parallel for
        for (size_t j = 0; j < sphereSet.total_combinations; ++j)
        {
            Scatterer scatterer = sphereSet.get_scatterer_by_index(j, source);

            #pragma omp parallel for
            for (size_t k = 0; k < detectorSet.total_combinations; ++k)
            {
                Detector detector = detectorSet.get_detector_by_index(k);
                size_t idx = flatten_multi_index(source.indices, scatterer.indices, detector.indices, array_shape);

                output_array[idx] = detector.get_coupling(scatterer);
            }
        }
    }

    return vector_to_numpy(output_array, array_shape);
}