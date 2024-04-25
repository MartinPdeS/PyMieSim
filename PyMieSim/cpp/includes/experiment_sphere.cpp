#pragma once

template<typename Function>
pybind11::array_t<complex128> Experiment::get_sphere_coefficient(Function function, size_t max_order) const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<complex128> output_array(full_size);

    #pragma omp parallel for collapse(5)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t sd=0; sd<array_shape[2]; ++sd)
    for (size_t si=0; si<array_shape[3]; ++si)
    for (size_t mi=0; mi<array_shape[4]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, sd, si, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        SPHERE::Scatterer scatterer = sphereSet.to_object(sd, si, wl, mi, source);

        output_array[idx] = (scatterer.*function)()[max_order];
    }

  return vector_to_numpy(output_array, array_shape);
}

template<typename Function>
pybind11::array_t<double> Experiment::get_sphere_data(Function function) const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(5)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t sd=0; sd<array_shape[2]; ++sd)
    for (size_t si=0; si<array_shape[3]; ++si)
    for (size_t mi=0; mi<array_shape[4]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, sd, si, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        SPHERE::Scatterer scatterer = sphereSet.to_object(sd, si, wl, mi, source);

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


pybind11::array_t<double> Experiment::get_sphere_coupling() const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(11)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t sd=0; sd<array_shape[2]; ++sd)
    for (size_t si=0; si<array_shape[3]; ++si)
    for (size_t mi=0; mi<array_shape[4]; ++mi)
    for (size_t sf=0; sf<array_shape[5]; ++sf)
    for (size_t ra=0; ra<array_shape[6]; ++ra)
    for (size_t na=0; na<array_shape[7]; ++na)
    for (size_t po=0; po<array_shape[8]; ++po)
    for (size_t go=0; go<array_shape[9]; ++go)
    for (size_t pf=0; pf<array_shape[10]; ++pf)
    {
        size_t idx = flatten_multi_index({wl, jv, sd, si, mi, sf, ra, na, po, go, pf}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        DETECTOR::Detector detector = detectorSet.to_object(sf, ra, na, po, go, pf);

        SPHERE::Scatterer scatterer = sphereSet.to_object(sd, si, wl, mi, source);

        double coupling = detector.get_coupling(scatterer);

        output_array[idx] = abs(coupling);
    }

    return vector_to_numpy(output_array, array_shape);
}
