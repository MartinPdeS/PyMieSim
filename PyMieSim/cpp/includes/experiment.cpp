#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "definitions.cpp"
#include "numpy_interface.cpp"
#include "sets.cpp"
#include <iostream>


class Experiment
{
    public:
        SPHERE::Set sphereSet;
        CYLINDER::Set cylinderSet;
        CORESHELL::Set coreshellSet;
        DETECTOR::Set detectorSet;
        SOURCE::Set sourceSet;

        Experiment(){}

        void set_sphere(SPHERE::Set& ScattererSet){this->sphereSet = ScattererSet;}

        void set_cylinder(CYLINDER::Set& ScattererSet) {this->cylinderSet = ScattererSet;}

        void set_coreshell(CORESHELL::Set& ScattererSet) {this->coreshellSet = ScattererSet;}

        void set_source(SOURCE::Set &sourceSet){this->sourceSet = sourceSet;}

        void set_detector(DETECTOR::Set &detectorSet){this->detectorSet = detectorSet;}


        //--------------------------------------SPHERE------------------------------------
        pybind11::array_t<complex128>
        get_sphere_coefficient(std::vector<complex128> (SPHERE::Scatterer::*function)(void), size_t max_order=0)
        {
            if (sphereSet.bounded_index)
                return get_sphere_coefficient_material(function, max_order);
            else
                return get_sphere_coefficient_index(function, max_order);
        }

        pybind11::array_t<double>
        get_sphere_data(double (SPHERE::Scatterer::*function)(void), size_t max_order=0)
        {
            if (sphereSet.bounded_index)
                return get_sphere_data_material(function, max_order);
            else
                return get_sphere_data_index(function, max_order);
        }


        pybind11::array_t<double>
        get_sphere_coupling()
        {
            if (sphereSet.bounded_index)
                return get_sphere_coupling_material();
            else
                return get_sphere_coupling_index();
        }


        size_t flatten_multi_index(const std::vector<size_t> &multi_index, const std::vector<size_t> &dimension)
        {
            size_t
                flatten_index = 0,
                multiplier = 1,
                dimension_size,
                index;

            for (size_t index_number=0; index_number<multi_index.size(); index_number++)
            {
                index = multi_index[index_number];
                multiplier = 1;
                for (size_t dim_number=index_number+1; dim_number<dimension.size(); dim_number++)
                {
                    dimension_size = dimension[dim_number];
                    multiplier *= dimension_size;
                }

                flatten_index += index * multiplier;
            }

            return flatten_index;
        }

        pybind11::array_t<complex128>
        get_sphere_coefficient_material(std::vector<complex128> (SPHERE::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace SPHERE;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                sphereSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<complex128> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                      sourceSet.wavelength[w],
                                      sourceSet.jones_vector[j],
                                      sourceSet.amplitude[w]
                                    );

                                    SPHERE::State scatterer_state = SPHERE::State(
                                      sphereSet.diameter[d],
                                      sphereSet.material[i][w],
                                      sphereSet.n_medium[n]
                                    );

                                    SPHERE::Scatterer scatterer = SPHERE::Scatterer(scatterer_state, source_state, max_order+1);

                                    output_array[idx] = 0. * (scatterer.*function)()[max_order];
                                }

          return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<complex128>
        get_sphere_coefficient_index(std::vector<complex128> (SPHERE::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace SPHERE;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                sphereSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<complex128> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                        sourceSet.wavelength[w],
                                        sourceSet.jones_vector[j],
                                        sourceSet.amplitude[w]
                                    );

                                    SPHERE::State scatterer_state = SPHERE::State(
                                        sphereSet.diameter[d],
                                        sphereSet.index[i],
                                        sphereSet.n_medium[n]
                                    );


                                    SPHERE::Scatterer scatterer = SPHERE::Scatterer(
                                        scatterer_state,
                                        source_state,
                                        max_order+1
                                    );

                                    output_array[idx] = (scatterer.*function)()[max_order];
                                }

          return vector_to_ndarray(output_array, array_shape);
        }






        pybind11::array_t<double>
        get_sphere_data_material(double (SPHERE::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace SPHERE;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                sphereSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                        sourceSet.wavelength[w],
                                        sourceSet.jones_vector[j],
                                        sourceSet.amplitude[w]
                                    );

                                    SPHERE::State scatterer_state = SPHERE::State(
                                        sphereSet.diameter[d],
                                        sphereSet.material[i][w],
                                        sphereSet.n_medium[n]
                                    );


                                    SPHERE::Scatterer scatterer = SPHERE::Scatterer(
                                        scatterer_state,
                                        source_state
                                    );

                                    output_array[idx] = (scatterer.*function)();
                                }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double>
        get_sphere_data_index(double (SPHERE::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace SPHERE;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                sphereSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(6) shared(output_array)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                        sourceSet.wavelength[w],
                                        sourceSet.jones_vector[j],
                                        sourceSet.amplitude[w]
                                    );

                                    SPHERE::State scatterer_state = SPHERE::State(
                                        sphereSet.diameter[d],
                                        sphereSet.index[i],
                                        sphereSet.n_medium[n]
                                    );


                                    SPHERE::Scatterer Scat = SPHERE::Scatterer(
                                        scatterer_state,
                                        source_state,
                                        max_order
                                    );

                                    output_array[idx] = (Scat.*function)();
                                }

            return vector_to_ndarray(output_array, array_shape);
        }






        pybind11::array_t<double> get_sphere_coupling_material()
        {
            using namespace SPHERE;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                sphereSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(11)
                for (size_t w=0; w<array_shape[0]; ++w)
                    for (size_t j=0; j<array_shape[1]; ++j)
                            for (size_t d=0; d<array_shape[2]; ++d)
                                for (size_t i=0; i<array_shape[3]; ++i)
                                    for (size_t n=0; n<array_shape[4]; ++n)
                                        for (size_t s=0; s<array_shape[5]; ++s)
                                            for (size_t na=0; na<array_shape[6]; ++na)
                                                for (size_t p=0; p<array_shape[7]; ++p)
                                                    for (size_t g=0; g<array_shape[8]; ++g)
                                                        for (size_t f=0; f<array_shape[9]; ++f)
                                                        {
                                                            size_t idx = flatten_multi_index({w, j, d, i, n, s, na, p, g, f}, array_shape);

                                                            SOURCE::State source_state = SOURCE::State(
                                                                sourceSet.wavelength[w],
                                                                sourceSet.jones_vector[j],
                                                                sourceSet.amplitude[w]
                                                            );

                                                            SPHERE::State scatterer_state = SPHERE::State(
                                                                sphereSet.diameter[d],
                                                                sphereSet.material[i][w],
                                                                sphereSet.n_medium[n]
                                                            );


                                                            DETECTOR::State detector_state  = DETECTOR::State(
                                                                detectorSet.scalar_field[s],
                                                                detectorSet.NA[na],
                                                                detectorSet.phi_offset[p],
                                                                detectorSet.gamma_offset[g],
                                                                detectorSet.polarization_filter[f],
                                                                detectorSet.coherent,
                                                                detectorSet.point_coupling
                                                            );

                                                            SPHERE::Scatterer scatterer = SPHERE::Scatterer(
                                                                scatterer_state,
                                                                source_state
                                                            );

                                                            DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                            output_array[idx] = abs( detector.Coupling(scatterer) );
                                                        }

            return vector_to_ndarray(output_array, array_shape);
        }



        pybind11::array_t<double> get_sphere_coupling_index()
        {
            using namespace SPHERE;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                sphereSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(11)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                    for (size_t s=0; s<array_shape[5]; ++s)
                                        for (size_t na=0; na<array_shape[6]; ++na)
                                            for (size_t p=0; p<array_shape[7]; ++p)
                                                for (size_t g=0; g<array_shape[8]; ++g)
                                                    for (size_t f=0; f<array_shape[9]; ++f)
                                                    {
                                                        size_t idx = flatten_multi_index({w, j, d, i, n, s, na, p, g, f}, array_shape);


                                                        SOURCE::State source_state = SOURCE::State(
                                                            sourceSet.wavelength[w],
                                                            sourceSet.jones_vector[j],
                                                            sourceSet.amplitude[w]
                                                        );

                                                        SPHERE::State scatterer_state = SPHERE::State(
                                                            sphereSet.diameter[d],
                                                            sphereSet.index[i],
                                                            sphereSet.n_medium[n]
                                                        );


                                                        DETECTOR::State detector_state = DETECTOR::State(
                                                            detectorSet.scalar_field[s],
                                                            detectorSet.NA[na],
                                                            detectorSet.phi_offset[p],
                                                            detectorSet.gamma_offset[g],
                                                            detectorSet.polarization_filter[f],
                                                            detectorSet.coherent,
                                                            detectorSet.point_coupling
                                                        );

                                                        SPHERE::Scatterer scatterer = SPHERE::Scatterer(
                                                            scatterer_state,
                                                            source_state
                                                        );

                                                        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                        double coupling = detector.Coupling(scatterer);

                                                        output_array[idx] = abs(coupling);
                                                    }

            return vector_to_ndarray(output_array, array_shape);
        }


        //--------------------------------------CYLINDER------------------------------------
        pybind11::array_t<complex128> get_cylinder_coefficient(std::vector<complex128> (CYLINDER::Scatterer::*function)(void), size_t max_order=0)
        {
            if (cylinderSet.bounded_index)
                return get_cylinder_coefficient_material(function, max_order);
            else
                return get_cylinder_coefficient_index(function, max_order);
        }

        pybind11::array_t<double> get_cylinder_data(double (CYLINDER::Scatterer::*function)(void), size_t max_order=0)
        {
            if (cylinderSet.bounded_index)
                return get_cylinder_data_material(function, max_order);
            else
                return get_cylinder_data_index(function, max_order);
        }

        pybind11::array_t<double> get_cylinder_coupling()
        {
            if (cylinderSet.bounded_index)
                return get_cylinder_coupling_bound();
            else
                return get_cylinder_coupling_unbound();
        }



        pybind11::array_t<complex128> get_cylinder_coefficient_material(std::vector<complex128> (CYLINDER::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CYLINDER;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                cylinderSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<complex128> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                        sourceSet.wavelength[w],
                                        sourceSet.jones_vector[j],
                                        sourceSet.amplitude[w]
                                    );

                                    CYLINDER::State scatterer_state = CYLINDER::State(
                                        cylinderSet.diameter[d],
                                        cylinderSet.material[i][w],
                                        cylinderSet.n_medium[n]
                                    );

                                    CYLINDER::Scatterer Scat = CYLINDER::Scatterer(
                                        scatterer_state,
                                        source_state,
                                        max_order + 1
                                    );

                                    output_array[idx] = (Scat.*function)()[max_order];
                                }

            return vector_to_ndarray(output_array, array_shape);
        }






        pybind11::array_t<complex128> get_cylinder_coefficient_index(std::vector<complex128> (CYLINDER::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CYLINDER;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                cylinderSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<complex128> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                        sourceSet.wavelength[w],
                                        sourceSet.jones_vector[j],
                                        sourceSet.amplitude[w]
                                    );

                                    CYLINDER::State scatterer_state = CYLINDER::State(
                                        cylinderSet.diameter[d],
                                        cylinderSet.index[i],
                                        cylinderSet.n_medium[n]
                                    );

                                    CYLINDER::Scatterer Scat = CYLINDER::Scatterer(
                                        scatterer_state,
                                        source_state,
                                        max_order + 1
                                    );

                                    output_array[idx] = (Scat.*function)()[max_order];
                                }

            return vector_to_ndarray(output_array, array_shape);
        }



        pybind11::array_t<double> get_cylinder_data_material(double (CYLINDER::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CYLINDER;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                cylinderSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                        sourceSet.wavelength[w],
                                        sourceSet.jones_vector[j],
                                        sourceSet.amplitude[w]
                                    );

                                    CYLINDER::State scatterer_state = CYLINDER::State(
                                        cylinderSet.diameter[d],
                                        cylinderSet.material[i][w],
                                        cylinderSet.n_medium[n]
                                    );

                                    CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
                                        scatterer_state,
                                        source_state,
                                        max_order
                                    );

                                    output_array[idx] = (scatterer.*function)();
                                }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double> get_cylinder_data_index(double (CYLINDER::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CYLINDER;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                cylinderSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);


            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                {
                                    size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

                                    SOURCE::State source_state = SOURCE::State(
                                        sourceSet.wavelength[w],
                                        sourceSet.jones_vector[j],
                                        sourceSet.amplitude[w]
                                    );

                                    CYLINDER::State scatterer_state = CYLINDER::State(
                                        cylinderSet.diameter[d],
                                        cylinderSet.index[i],
                                        cylinderSet.n_medium[n]
                                    );

                                    CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
                                        scatterer_state,
                                        source_state,
                                        max_order
                                    );

                                    output_array[idx] = (scatterer.*function)();
                                }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double> get_cylinder_coupling_bound()
        {
            using namespace CYLINDER;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                cylinderSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(11)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                    for (size_t s=0; s<array_shape[5]; ++s)
                                        for (size_t na=0; na<array_shape[6]; ++na)
                                            for (size_t p=0; p<array_shape[7]; ++p)
                                                for (size_t g=0; g<array_shape[8]; ++g)
                                                    for (size_t f=0; f<array_shape[9]; ++f)
                                                    {
                                                        size_t idx = flatten_multi_index({w, j, d, i, n, s, na, p, g, f}, array_shape);


                                                        SOURCE::State source_state = SOURCE::State(
                                                            sourceSet.wavelength[w],
                                                            sourceSet.jones_vector[j],
                                                            sourceSet.amplitude[w]
                                                        );

                                                        CYLINDER::State scatterer_state = CYLINDER::State(
                                                            cylinderSet.diameter[d],
                                                            cylinderSet.material[i][w],
                                                            cylinderSet.n_medium[n]
                                                        );

                                                        DETECTOR::State detector_state = DETECTOR::State(
                                                            detectorSet.scalar_field[s],
                                                            detectorSet.NA[na],
                                                            detectorSet.phi_offset[p],
                                                            detectorSet.gamma_offset[g],
                                                            detectorSet.polarization_filter[f],
                                                            detectorSet.coherent,
                                                            detectorSet.point_coupling
                                                        );

                                                        CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
                                                            scatterer_state,
                                                            source_state
                                                        );

                                                        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                        output_array[idx] = abs( detector.Coupling(scatterer) );
                                                    }

            return vector_to_ndarray(output_array, array_shape);
        }






        pybind11::array_t<double> get_cylinder_coupling_unbound()
        {
            using namespace CYLINDER;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                cylinderSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(11)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t d=0; d<array_shape[2]; ++d)
                            for (size_t i=0; i<array_shape[3]; ++i)
                                for (size_t n=0; n<array_shape[4]; ++n)
                                    for (size_t s=0; s<array_shape[5]; ++s)
                                        for (size_t na=0; na<array_shape[6]; ++na)
                                            for (size_t p=0; p<array_shape[7]; ++p)
                                                for (size_t g=0; g<array_shape[8]; ++g)
                                                    for (size_t f=0; f<array_shape[9]; ++f)
                                                    {
                                                        size_t idx = flatten_multi_index({w, j, d, i, n, s, na, p, g, f}, array_shape);

                                                        SOURCE::State source_state = SOURCE::State(
                                                            sourceSet.wavelength[w],
                                                            sourceSet.jones_vector[j],
                                                            sourceSet.amplitude[w]
                                                        );

                                                        CYLINDER::State scatterer_state = CYLINDER::State(
                                                            cylinderSet.diameter[d],
                                                            cylinderSet.index[i],
                                                            cylinderSet.n_medium[n]
                                                        );

                                                        DETECTOR::State detector_state = DETECTOR::State(
                                                            detectorSet.scalar_field[s],
                                                            detectorSet.NA[na],
                                                            detectorSet.phi_offset[p],
                                                            detectorSet.gamma_offset[g],
                                                            detectorSet.polarization_filter[f],
                                                            detectorSet.coherent,
                                                            detectorSet.point_coupling
                                                        );

                                                        CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
                                                            scatterer_state,
                                                            source_state
                                                        );

                                                        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                        output_array[idx] = abs(detector.Coupling(scatterer));
                                                    }

            return vector_to_ndarray(output_array, array_shape);
        }




        //--------------------------------------CORESHELL------------------------------------
        pybind11::array_t<complex128> get_coreshell_coefficient(std::vector<complex128> (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            if (coreshellSet.bounded_core && coreshellSet.bounded_shell)
                return get_coreshell_coefficient_core_material_shell_material(function, max_order);

            if (coreshellSet.bounded_core && !coreshellSet.bounded_shell)
                return get_coreshell_coefficient_core_material_shell_index(function, max_order);

            if (!coreshellSet.bounded_core && coreshellSet.bounded_shell)
                return get_coreshell_coefficient_core_index_shell_material(function, max_order);

            if (!coreshellSet.bounded_core && !coreshellSet.bounded_shell)
                return get_coreshell_coefficient_core_index_shell_index(function, max_order);
        }

        pybind11::array_t<double> get_coreshell_data(double (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            if (coreshellSet.bounded_core && coreshellSet.bounded_shell)
                return get_coreshell_data_core_material_shell_material(function, max_order);

            if (coreshellSet.bounded_core && !coreshellSet.bounded_shell)
                return get_coreshell_data_core_material_shell_index(function, max_order);

            if (!coreshellSet.bounded_core && coreshellSet.bounded_shell)
                return get_coreshell_data_core_index_shell_material(function, max_order);

            if (!coreshellSet.bounded_core && !coreshellSet.bounded_shell)
                return get_coreshell_data_core_index_shell_index(function, max_order);
        }


        pybind11::array_t<double> get_coreshell_coupling()
        {
            if (coreshellSet.bounded_core && coreshellSet.bounded_shell)
                return get_coreshell_coupling_core_material_shell_material();

            if (coreshellSet.bounded_core && !coreshellSet.bounded_shell)
                return get_coreshell_coupling_core_material_shell_index();

            if (!coreshellSet.bounded_core && coreshellSet.bounded_shell)
                return get_coreshell_coupling_core_index_shell_material();

            if (!coreshellSet.bounded_core && !coreshellSet.bounded_shell)
                return get_coreshell_coupling_core_index_shell_index();
        }


        pybind11::array_t<complex128> get_coreshell_coefficient_core_material_shell_index(std::vector<complex128> (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<complex128> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                            for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                                for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                    for (size_t Si=0; Si<array_shape[5]; ++Si)
                                        for (size_t n=0; n<array_shape[6]; ++n)
                                        {
                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

                                            SOURCE::State source_state = SOURCE::State(
                                                sourceSet.wavelength[w],
                                                sourceSet.jones_vector[j],
                                                sourceSet.amplitude[w]
                                            );

                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                coreshellSet.core_diameter[Cd],
                                                coreshellSet.shell_width[Sd],
                                                coreshellSet.core_material[Ci][w],
                                                coreshellSet.shell_index[Si],
                                                coreshellSet.n_medium[n]
                                            );

                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                scatterer_state,
                                                source_state,
                                                max_order + 1
                                            );

                                            output_array[idx] = (scatterer.*function)()[max_order];
                                        }

            return vector_to_ndarray(output_array, array_shape);
        }



        pybind11::array_t<complex128> get_coreshell_coefficient_core_index_shell_material(std::vector<complex128> (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<complex128> output_array(full_size);

            #pragma omp parallel for collapse(8)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                            for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                                for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                    for (size_t Si=0; Si<array_shape[5]; ++Si)
                                        for (size_t n=0; n<array_shape[6]; ++n)
                                        {
                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

                                            SOURCE::State source_state = SOURCE::State(
                                                sourceSet.wavelength[w],
                                                sourceSet.jones_vector[j],
                                                sourceSet.amplitude[w]
                                            );

                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                coreshellSet.core_diameter[Cd],
                                                coreshellSet.shell_width[Sd],
                                                coreshellSet.core_index[Ci],
                                                coreshellSet.shell_material[Si][w],
                                                coreshellSet.n_medium[n]
                                            );

                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                scatterer_state,
                                                source_state,
                                                max_order + 1
                                            );

                                            output_array[idx] = (scatterer.*function)()[max_order];
                                        }

            return vector_to_ndarray(output_array, array_shape);
        }



        pybind11::array_t<complex128> get_coreshell_coefficient_core_material_shell_material(std::vector<complex128> (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
        using namespace CORESHELL;

        std::vector<size_t> array_shape = concatenate_vector(
            sourceSet.get_array_shape(),
            coreshellSet.get_array_shape()
        );

        size_t full_size = get_vector_sigma(array_shape);

        std::vector<complex128> output_array(full_size);

        #pragma omp parallel for collapse(8)
        for (size_t w=0; w<array_shape[0]; ++w)
            for (size_t j=0; j<array_shape[1]; ++j)
                    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                        for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                            for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                for (size_t Si=0; Si<array_shape[5]; ++Si)
                                    for (size_t n=0; n<array_shape[6]; ++n)
                                    {
                                        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

                                        SOURCE::State source_state = SOURCE::State(
                                            sourceSet.wavelength[w],
                                            sourceSet.jones_vector[j],
                                            sourceSet.amplitude[w]
                                        );

                                        CORESHELL::State scatterer_state = CORESHELL::State(
                                            coreshellSet.core_diameter[Cd],
                                            coreshellSet.shell_width[Sd],
                                            coreshellSet.core_material[Ci][w],
                                            coreshellSet.shell_material[Si][w],
                                            coreshellSet.n_medium[n]
                                        );

                                        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                            scatterer_state,
                                            source_state,
                                            max_order + 1
                                        );

                                        output_array[idx] = (scatterer.*function)()[max_order];
                                    }

        return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<complex128> get_coreshell_coefficient_core_index_shell_index(std::vector<complex128> (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<complex128> output_array(full_size);

            #pragma omp parallel for collapse(8)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                            for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                                for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                    for (size_t Si=0; Si<array_shape[5]; ++Si)
                                        for (size_t n=0; n<array_shape[6]; ++n)
                                        {
                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

                                            SOURCE::State source_state = SOURCE::State(
                                                sourceSet.wavelength[w],
                                                sourceSet.jones_vector[j],
                                                sourceSet.amplitude[w]
                                            );

                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                coreshellSet.core_diameter[Cd],
                                                coreshellSet.shell_width[Sd],
                                                coreshellSet.core_index[Ci],
                                                coreshellSet.shell_index[Si],
                                                coreshellSet.n_medium[n]
                                            );

                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                scatterer_state,
                                                source_state,
                                                max_order + 1
                                            );

                                            output_array[idx] = (scatterer.*function)()[max_order];
                                        }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double> get_coreshell_data_core_material_shell_index(double (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(6)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                            for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                                for (size_t Cm=0; Cm<array_shape[4]; ++Cm)
                                    for (size_t Si=0; Si<array_shape[5]; ++Si)
                                        for (size_t n=0; n<array_shape[6]; ++n)
                                        {
                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Cm, Si, n}, array_shape);

                                            SOURCE::State source_state = SOURCE::State(
                                                sourceSet.wavelength[w],
                                                sourceSet.jones_vector[j],
                                                sourceSet.amplitude[w]
                                            );

                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                coreshellSet.core_diameter[Cd],
                                                coreshellSet.shell_width[Sd],
                                                coreshellSet.core_material[Cm][w],
                                                coreshellSet.shell_index[Si],
                                                coreshellSet.n_medium[n]
                                            );

                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                scatterer_state,
                                                source_state,
                                                max_order
                                            );

                                            output_array[idx] = (scatterer.*function)();
                                        }

            return vector_to_ndarray(output_array, array_shape);
        }



        pybind11::array_t<double> get_coreshell_data_core_index_shell_material(double (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(8)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                        for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                            for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                for (size_t Sm=0; Sm<array_shape[5]; ++Sm)
                                    for (size_t n=0; n<array_shape[6]; ++n)
                                    {
                                        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Sm, n}, array_shape);

                                        SOURCE::State source_state = SOURCE::State(
                                            sourceSet.wavelength[w],
                                            sourceSet.jones_vector[j],
                                            sourceSet.amplitude[w]
                                        );

                                        CORESHELL::State scatterer_state = CORESHELL::State(
                                            coreshellSet.core_diameter[Cd],
                                            coreshellSet.shell_width[Sd],
                                            coreshellSet.core_index[Ci],
                                            coreshellSet.shell_material[Sm][w],
                                            coreshellSet.n_medium[n]
                                        );

                                        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                            scatterer_state,
                                            source_state,
                                            max_order
                                        );

                                        output_array[idx] = (scatterer.*function)();
                                    }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double> get_coreshell_data_core_material_shell_material(double (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(8)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                        for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                            for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                                for (size_t Cm=0; Cm<array_shape[4]; ++Cm)
                                    for (size_t Sm=0; Sm<array_shape[5]; ++Sm)
                                        for (size_t n=0; n<array_shape[6]; ++n)
                                        {
                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Cm, Sm, n}, array_shape);

                                            SOURCE::State source_state = SOURCE::State(
                                                sourceSet.wavelength[w],
                                                sourceSet.jones_vector[j],
                                                sourceSet.amplitude[w]
                                            );

                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                coreshellSet.core_diameter[Cd],
                                                coreshellSet.shell_width[Sd],
                                                coreshellSet.core_material[Cm][w],
                                                coreshellSet.shell_material[Sm][w],
                                                coreshellSet.n_medium[n]
                                            );

                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                scatterer_state,
                                                source_state,
                                                max_order
                                            );

                                            output_array[idx] = (scatterer.*function)();
                                        }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double> get_coreshell_data_core_index_shell_index(double (CORESHELL::Scatterer::*function)(void), size_t max_order=0)
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(8)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                        for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                            for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                for (size_t Si=0; Si<array_shape[5]; ++Si)
                                    for (size_t n=0; n<array_shape[6]; ++n)
                                    {
                                        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

                                        SOURCE::State source_state = SOURCE::State(
                                            sourceSet.wavelength[w],
                                            sourceSet.jones_vector[j],
                                            sourceSet.amplitude[w]
                                        );

                                        CORESHELL::State scatterer_state = CORESHELL::State(
                                            coreshellSet.core_diameter[Cd],
                                            coreshellSet.shell_width[Sd],
                                            coreshellSet.core_index[Ci],
                                            coreshellSet.shell_index[Si],
                                            coreshellSet.n_medium[n]
                                        );

                                        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                            scatterer_state,
                                            source_state,
                                            max_order
                                        );

                                        output_array[idx] = (scatterer.*function)();
                                    }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double> get_coreshell_coupling_core_index_shell_index()
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);


            #pragma omp parallel for collapse(13)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                        for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                            for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                for (size_t Si=0; Si<array_shape[5]; ++Si)
                                    for (size_t n=0; n<array_shape[6]; ++n)
                                        for (size_t s=0; s<array_shape[7]; ++s)
                                            for (size_t na=0; na<array_shape[8]; ++na)
                                                for (size_t p=0; p<array_shape[9]; ++p)
                                                    for (size_t g=0; g<array_shape[10]; ++g)
                                                        for (size_t f=0; f<array_shape[11]; ++f)
                                                        {
                                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

                                                            SOURCE::State source_state = SOURCE::State(
                                                                sourceSet.wavelength[w],
                                                                sourceSet.jones_vector[j],
                                                                sourceSet.amplitude[w]
                                                            );

                                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                                coreshellSet.core_diameter[Cd],
                                                                coreshellSet.shell_width[Sd],
                                                                coreshellSet.core_index[Ci],
                                                                coreshellSet.shell_index[Si],
                                                                coreshellSet.n_medium[n]
                                                            );

                                                            DETECTOR::State detector_state = DETECTOR::State(
                                                                detectorSet.scalar_field[s],
                                                                detectorSet.NA[na],
                                                                detectorSet.phi_offset[p],
                                                                detectorSet.gamma_offset[g],
                                                                detectorSet.polarization_filter[f],
                                                                detectorSet.coherent,
                                                                detectorSet.point_coupling
                                                            );

                                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                                scatterer_state,
                                                                source_state
                                                            );

                                                            DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                            output_array[idx] = abs( detector.Coupling(scatterer) );
                                                        }

            return vector_to_ndarray(output_array, array_shape);
        }




        pybind11::array_t<double> get_coreshell_coupling_core_material_shell_index()
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(13)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                        for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                            for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                for (size_t Si=0; Si<array_shape[5]; ++Si)
                                    for (size_t n=0; n<array_shape[6]; ++n)
                                        for (size_t s=0; s<array_shape[7]; ++s)
                                            for (size_t na=0; na<array_shape[8]; ++na)
                                                for (size_t p=0; p<array_shape[9]; ++p)
                                                    for (size_t g=0; g<array_shape[10]; ++g)
                                                        for (size_t f=0; f<array_shape[11]; ++f)
                                                        {
                                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

                                                            SOURCE::State source_state = SOURCE::State(
                                                                sourceSet.wavelength[w],
                                                                sourceSet.jones_vector[j],
                                                                sourceSet.amplitude[w]
                                                            );

                                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                                coreshellSet.core_diameter[Cd],
                                                                coreshellSet.shell_width[Sd],
                                                                coreshellSet.core_material[Ci][w],
                                                                coreshellSet.shell_index[Si],
                                                                coreshellSet.n_medium[n]
                                                            );

                                                            DETECTOR::State detector_state = DETECTOR::State(
                                                                detectorSet.scalar_field[s],
                                                                detectorSet.NA[na],
                                                                detectorSet.phi_offset[p],
                                                                detectorSet.gamma_offset[g],
                                                                detectorSet.polarization_filter[f],
                                                                detectorSet.coherent,
                                                                detectorSet.point_coupling
                                                            );

                                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                                scatterer_state,
                                                                source_state
                                                            );

                                                            DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                            output_array[idx] = abs( detector.Coupling(scatterer) );
                                                        }

            return vector_to_ndarray(output_array, array_shape);
        }




        pybind11::array_t<double> get_coreshell_coupling_core_index_shell_material()
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(13)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                        for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                            for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                for (size_t Si=0; Si<array_shape[5]; ++Si)
                                    for (size_t n=0; n<array_shape[6]; ++n)
                                        for (size_t s=0; s<array_shape[7]; ++s)
                                            for (size_t na=0; na<array_shape[8]; ++na)
                                                for (size_t p=0; p<array_shape[9]; ++p)
                                                    for (size_t g=0; g<array_shape[10]; ++g)
                                                        for (size_t f=0; f<array_shape[11]; ++f)
                                                        {
                                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

                                                            SOURCE::State source_state = SOURCE::State(
                                                                sourceSet.wavelength[w],
                                                                sourceSet.jones_vector[j],
                                                                sourceSet.amplitude[w]
                                                            );

                                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                                coreshellSet.core_diameter[Cd],
                                                                coreshellSet.shell_width[Sd],
                                                                coreshellSet.core_index[Ci],
                                                                coreshellSet.shell_material[Si][w],
                                                                coreshellSet.n_medium[n]
                                                            );

                                                            DETECTOR::State detector_state = DETECTOR::State(
                                                                detectorSet.scalar_field[s],
                                                                detectorSet.NA[na],
                                                                detectorSet.phi_offset[p],
                                                                detectorSet.gamma_offset[g],
                                                                detectorSet.polarization_filter[f],
                                                                detectorSet.coherent,
                                                                detectorSet.point_coupling
                                                            );

                                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                                scatterer_state,
                                                                source_state
                                                            );

                                                            DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                            output_array[idx] = abs( detector.Coupling(scatterer) );
                                                        }

            return vector_to_ndarray(output_array, array_shape);
        }


        pybind11::array_t<double> get_coreshell_coupling_core_material_shell_material()
        {
            using namespace CORESHELL;

            std::vector<size_t> array_shape = concatenate_vector(
                sourceSet.get_array_shape(),
                coreshellSet.get_array_shape(),
                detectorSet.get_array_shape()
            );

            size_t full_size = get_vector_sigma(array_shape);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for collapse(13)
            for (size_t w=0; w<array_shape[0]; ++w)
                for (size_t j=0; j<array_shape[1]; ++j)
                    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
                        for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
                            for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
                                for (size_t Si=0; Si<array_shape[5]; ++Si)
                                    for (size_t n=0; n<array_shape[6]; ++n)
                                        for (size_t s=0; s<array_shape[7]; ++s)
                                            for (size_t na=0; na<array_shape[8]; ++na)
                                                for (size_t p=0; p<array_shape[9]; ++p)
                                                    for (size_t g=0; g<array_shape[10]; ++g)
                                                        for (size_t f=0; f<array_shape[11]; ++f)
                                                        {
                                                            size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

                                                            SOURCE::State source_state = SOURCE::State(
                                                                sourceSet.wavelength[w],
                                                                sourceSet.jones_vector[j],
                                                                sourceSet.amplitude[w]
                                                            );

                                                            CORESHELL::State scatterer_state = CORESHELL::State(
                                                                coreshellSet.core_diameter[Cd],
                                                                coreshellSet.shell_width[Sd],
                                                                coreshellSet.core_material[Ci][w],
                                                                coreshellSet.shell_material[Si][w],
                                                                coreshellSet.n_medium[n]
                                                            );

                                                            DETECTOR::State  detector_state = DETECTOR::State(
                                                                detectorSet.scalar_field[s],
                                                                detectorSet.NA[na],
                                                                detectorSet.phi_offset[p],
                                                                detectorSet.gamma_offset[g],
                                                                detectorSet.polarization_filter[f],
                                                                detectorSet.coherent,
                                                                detectorSet.point_coupling
                                                            );

                                                            CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
                                                                scatterer_state,
                                                                source_state
                                                            );

                                                            DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

                                                            output_array[idx] = abs( detector.Coupling(scatterer) );
                                                        }

            return vector_to_ndarray(output_array, array_shape);
        }

pybind11::array_t<double> get_sphere_Qsca()    { return get_sphere_data( &SPHERE::Scatterer::get_Qsca ) ; }
pybind11::array_t<double> get_sphere_Qext()    { return get_sphere_data( &SPHERE::Scatterer::get_Qext ) ; }
pybind11::array_t<double> get_sphere_Qabs()    { return get_sphere_data( &SPHERE::Scatterer::get_Qabs ) ; }
pybind11::array_t<double> get_sphere_Qpr()     { return get_sphere_data( &SPHERE::Scatterer::get_Qpr ) ; }
pybind11::array_t<double> get_sphere_Qback()   { return get_sphere_data( &SPHERE::Scatterer::get_Qback ) ; }
pybind11::array_t<double> get_sphere_Qforward(){ return get_sphere_data( &SPHERE::Scatterer::get_Qforward ) ; }
pybind11::array_t<double> get_sphere_Csca()    { return get_sphere_data( &SPHERE::Scatterer::get_Csca ) ; }
pybind11::array_t<double> get_sphere_Cext()    { return get_sphere_data( &SPHERE::Scatterer::get_Cext ) ; }
pybind11::array_t<double> get_sphere_Cabs()    { return get_sphere_data( &SPHERE::Scatterer::get_Cabs ) ; }
pybind11::array_t<double> get_sphere_Cpr()     { return get_sphere_data( &SPHERE::Scatterer::get_Cpr ) ; }
pybind11::array_t<double> get_sphere_Cback()   { return get_sphere_data( &SPHERE::Scatterer::get_Cback ) ; }
pybind11::array_t<double> get_sphere_Cforward(){ return get_sphere_data( &SPHERE::Scatterer::get_Cforward ) ; }
pybind11::array_t<double> get_sphere_g()       { return get_sphere_data( &SPHERE::Scatterer::get_g ) ; }

pybind11::array_t<complex128> get_sphere_an(size_t max_order){ return get_sphere_coefficient( &SPHERE::Scatterer::get_an, max_order ) ; }
pybind11::array_t<complex128> get_sphere_bn(size_t max_order){ return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, max_order ) ; }
pybind11::array_t<complex128> get_sphere_a1()                { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 1 ) ; }
pybind11::array_t<complex128> get_sphere_b1()                { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 1 ) ; }
pybind11::array_t<complex128> get_sphere_a2()                { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 2 ) ; }
pybind11::array_t<complex128> get_sphere_b2()                { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 2 ) ; }
pybind11::array_t<complex128> get_sphere_a3()                { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 3 ) ; }
pybind11::array_t<complex128> get_sphere_b3()                { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 3 ) ; }

pybind11::array_t<double> get_cylinder_Qsca()    { return get_cylinder_data( &CYLINDER::Scatterer::get_Qsca ) ; }
pybind11::array_t<double> get_cylinder_Qext()    { return get_cylinder_data( &CYLINDER::Scatterer::get_Qext ) ; }
pybind11::array_t<double> get_cylinder_Qabs()    { return get_cylinder_data( &CYLINDER::Scatterer::get_Qabs ) ; }
pybind11::array_t<double> get_cylinder_Qpr()     { return get_cylinder_data( &CYLINDER::Scatterer::get_Qpr ) ; }
pybind11::array_t<double> get_cylinder_Qback()   { return get_cylinder_data( &CYLINDER::Scatterer::get_Qback ) ; }
pybind11::array_t<double> get_cylinder_Qforward(){ return get_cylinder_data( &CYLINDER::Scatterer::get_Qforward ) ; }
pybind11::array_t<double> get_cylinder_Csca()    { return get_cylinder_data( &CYLINDER::Scatterer::get_Csca ) ; }
pybind11::array_t<double> get_cylinder_Cext()    { return get_cylinder_data( &CYLINDER::Scatterer::get_Cext ) ; }
pybind11::array_t<double> get_cylinder_Cabs()    { return get_cylinder_data( &CYLINDER::Scatterer::get_Cabs ) ; }
pybind11::array_t<double> get_cylinder_Cpr()     { return get_cylinder_data( &CYLINDER::Scatterer::get_Cpr ) ; }
pybind11::array_t<double> get_cylinder_Cback()   { return get_cylinder_data( &CYLINDER::Scatterer::get_Cback ) ; }
pybind11::array_t<double> get_cylinder_Cforward(){ return get_cylinder_data( &CYLINDER::Scatterer::get_Cforward ) ; }
pybind11::array_t<double> get_cylinder_g()       { return get_cylinder_data( &CYLINDER::Scatterer::get_g ) ; }

pybind11::array_t<complex128> get_cylinder_a1n(size_t max_order){ return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, max_order ) ; }
pybind11::array_t<complex128> get_cylinder_b1n(size_t max_order){ return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, max_order ) ; }
pybind11::array_t<complex128> get_cylinder_a2n(size_t max_order){ return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, max_order ) ; }
pybind11::array_t<complex128> get_cylinder_b2n(size_t max_order){ return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, max_order ) ; }
pybind11::array_t<complex128> get_cylinder_a11()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, 1 ) ; }
pybind11::array_t<complex128> get_cylinder_b11()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, 1 ) ; }
pybind11::array_t<complex128> get_cylinder_a21()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, 1 ) ; }
pybind11::array_t<complex128> get_cylinder_b21()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, 1 ) ; }
pybind11::array_t<complex128> get_cylinder_a12()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, 2 ) ; }
pybind11::array_t<complex128> get_cylinder_b12()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, 2 ) ; }
pybind11::array_t<complex128> get_cylinder_a22()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, 2 ) ; }
pybind11::array_t<complex128> get_cylinder_b22()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, 2 ) ; }
pybind11::array_t<complex128> get_cylinder_a13()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, 3 ) ; }
pybind11::array_t<complex128> get_cylinder_b13()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, 3 ) ; }
pybind11::array_t<complex128> get_cylinder_a23()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, 3 ) ; }
pybind11::array_t<complex128> get_cylinder_b23()                { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, 3 ) ; }

pybind11::array_t<double> get_coreshell_Qsca()    { return get_coreshell_data( &CORESHELL::Scatterer::get_Qsca ) ; }
pybind11::array_t<double> get_coreshell_Qext()    { return get_coreshell_data( &CORESHELL::Scatterer::get_Qext ) ; }
pybind11::array_t<double> get_coreshell_Qabs()    { return get_coreshell_data( &CORESHELL::Scatterer::get_Qabs ) ; }
pybind11::array_t<double> get_coreshell_Qpr()     { return get_coreshell_data( &CORESHELL::Scatterer::get_Qpr ) ; }
pybind11::array_t<double> get_coreshell_Qback()   { return get_coreshell_data( &CORESHELL::Scatterer::get_Qback ) ; }
pybind11::array_t<double> get_coreshell_Qforward(){ return get_coreshell_data( &CORESHELL::Scatterer::get_Qforward ) ; }
pybind11::array_t<double> get_coreshell_Csca()    { return get_coreshell_data( &CORESHELL::Scatterer::get_Csca ) ; }
pybind11::array_t<double> get_coreshell_Cext()    { return get_coreshell_data( &CORESHELL::Scatterer::get_Cext ) ; }
pybind11::array_t<double> get_coreshell_Cabs()    { return get_coreshell_data( &CORESHELL::Scatterer::get_Cabs ) ; }
pybind11::array_t<double> get_coreshell_Cpr()     { return get_coreshell_data( &CORESHELL::Scatterer::get_Cpr ) ; }
pybind11::array_t<double> get_coreshell_Cback()   { return get_coreshell_data( &CORESHELL::Scatterer::get_Cback ) ; }
pybind11::array_t<double> get_coreshell_Cforward(){ return get_coreshell_data( &CORESHELL::Scatterer::get_Cforward ) ; }
pybind11::array_t<double> get_coreshell_g()       { return get_coreshell_data( &CORESHELL::Scatterer::get_g ) ; }

pybind11::array_t<complex128> get_coreshell_an(size_t max_order){ return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, max_order ) ; }
pybind11::array_t<complex128> get_coreshell_bn(size_t max_order){ return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, max_order ) ; }
pybind11::array_t<complex128> get_coreshell_a1()                { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 1 ) ; }
pybind11::array_t<complex128> get_coreshell_b1()                { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 1 ) ; }
pybind11::array_t<complex128> get_coreshell_a2()                { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 2 ) ; }
pybind11::array_t<complex128> get_coreshell_b2()                { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 2 ) ; }
pybind11::array_t<complex128> get_coreshell_a3()                { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 3 ) ; }
pybind11::array_t<complex128> get_coreshell_b3()                { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 3 ) ; }
};



#endif

