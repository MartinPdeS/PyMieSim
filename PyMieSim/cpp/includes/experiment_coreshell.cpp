#pragma once

template<typename Function>
pybind11::array_t<complex128> Experiment::get_coreshell_coefficient(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<complex128> output_array(full_size);

    #pragma omp parallel for collapse(9)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t na=0; na<array_shape[2]; ++na)
    for (size_t op=0; op<array_shape[3]; ++op)
    for (size_t cd=0; cd<array_shape[4]; ++cd)
    for (size_t sw=0; sw<array_shape[5]; ++sw)
    for (size_t cm=0; cm<array_shape[6]; ++cm)
    for (size_t sm=0; sm<array_shape[7]; ++sm)
    for (size_t mi=0; mi<array_shape[8]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, na, op, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Gaussian source = sourceSet.to_object(wl, jv, na, op);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

        output_array[idx] = (scatterer.*function)()[max_order];
    }

    return vector_to_numpy(output_array, array_shape);
}

template<typename Function>
pybind11::array_t<double> Experiment::get_coreshell_data(Function function) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(9)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t na=0; na<array_shape[2]; ++na)
    for (size_t op=0; op<array_shape[3]; ++op)
    for (size_t cd=0; cd<array_shape[4]; ++cd)
    for (size_t sw=0; sw<array_shape[5]; ++sw)
    for (size_t cm=0; cm<array_shape[6]; ++cm)
    for (size_t sm=0; sm<array_shape[7]; ++sm)
    for (size_t mi=0; mi<array_shape[8]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, na, op, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Gaussian source = sourceSet.to_object(wl, jv, na, op);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


pybind11::array_t<double> Experiment::get_coreshell_coupling() const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(16)
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t na_=0; na_<array_shape[2]; ++na_)
    for (size_t op=0; op<array_shape[3]; ++op)
    for (size_t cd=0; cd<array_shape[4]; ++cd)
    for (size_t sw=0; sw<array_shape[5]; ++sw)
    for (size_t cm=0; cm<array_shape[6]; ++cm)
    for (size_t sm=0; sm<array_shape[7]; ++sm)
    for (size_t mi=0; mi<array_shape[8]; ++mi)
    for (size_t mn=0; mn<array_shape[9]; ++mn)
    for (size_t fs=0; fs<array_shape[10]; ++fs)
    for (size_t ra=0; ra<array_shape[11]; ++ra)
    for (size_t na=0; na<array_shape[12]; ++na)
    for (size_t po=0; po<array_shape[13]; ++po)
    for (size_t go=0; go<array_shape[14]; ++go)
    for (size_t pf=0; pf<array_shape[15]; ++pf)
    {
        size_t idx = flatten_multi_index({wl, jv, na_, op, cd, sw, cm, sm, mi, mn, fs, ra, na, po, go, pf}, array_shape);

        SOURCE::Gaussian source = sourceSet.to_object(wl, jv, na_, op);

        DETECTOR::Detector detector = detectorSet.to_object(mn, fs, ra, na, po, go, pf);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}


