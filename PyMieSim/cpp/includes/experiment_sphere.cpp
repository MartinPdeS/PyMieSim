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

    #pragma omp parallel for collapse(7)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t na=0; na<array_shape[2]; ++na)
    for (size_t op=0; op<array_shape[3]; ++op)
    for (size_t sd=0; sd<array_shape[4]; ++sd)
    for (size_t si=0; si<array_shape[5]; ++si)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, na, op, sd, si, mi}, array_shape);

        SOURCE::Gaussian source = sourceSet.to_object(wl, jv, na, op);

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

    #pragma omp parallel for collapse(7)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t na=0; na<array_shape[2]; ++na)
    for (size_t op=0; op<array_shape[3]; ++op)
    for (size_t sd=0; sd<array_shape[4]; ++sd)
    for (size_t si=0; si<array_shape[5]; ++si)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, na, op, sd, si, mi}, array_shape);

        SOURCE::Gaussian source = sourceSet.to_object(wl, jv, na, op);

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

    #pragma omp parallel for collapse(14)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t na_=0; na_<array_shape[2]; ++na_)
    for (size_t op=0; op<array_shape[3]; ++op)
    for (size_t sd=0; sd<array_shape[4]; ++sd)
    for (size_t si=0; si<array_shape[5]; ++si)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    for (size_t mn=0; mn<array_shape[7]; ++mn)
    for (size_t fs=0; fs<array_shape[8]; ++fs)
    for (size_t ra=0; ra<array_shape[9]; ++ra)
    for (size_t na=0; na<array_shape[10]; ++na)
    for (size_t po=0; po<array_shape[11]; ++po)
    for (size_t go=0; go<array_shape[12]; ++go)
    for (size_t pf=0; pf<array_shape[13]; ++pf)
    {
        size_t idx = flatten_multi_index({wl, jv, na_, op, sd, si, mi, mn, fs, ra, na, po, go, pf}, array_shape);

        SOURCE::Gaussian source = sourceSet.to_object(wl, jv, na_, op);

        DETECTOR::Detector detector = detectorSet.to_object(mn, fs, ra, na, po, go, pf);

        SPHERE::Scatterer scatterer = sphereSet.to_object(sd, si, wl, mi, source);

        double coupling = detector.get_coupling(scatterer);

        output_array[idx] = abs(coupling);

    }

    return vector_to_numpy(output_array, array_shape);
}
