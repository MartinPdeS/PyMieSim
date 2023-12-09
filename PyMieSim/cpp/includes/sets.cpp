#ifndef SETS_H
#define SETS_H

#include "definitions.cpp"
#include "sources.cpp"
#include "sphere.cpp"
#include "cylinder.cpp"
#include "core_shell.cpp"
#include "detectors.cpp"

namespace SPHERE
{
    class Set
    {
        public:
            DVector diameter;
            CVector index;
            std::vector<std::vector<complex128>> material;
            DVector n_medium;
            bool bounded_index;
            std::vector<State> States;

            Set(){}
            Set(DVector &diameter, std::vector<std::vector<complex128>> &material, DVector &n_medium)
            : diameter(diameter), material(material), n_medium(n_medium)
            {
            bounded_index=true;
            }

            Set(DVector &diameter, std::vector<complex128> &index, DVector &n_medium)
            : diameter(diameter), index(index), n_medium(n_medium)
            {
            bounded_index=false;
            }

            State operator[](size_t idx){return this->States[idx];}

            std::vector<size_t> get_array_shape()
            {
                if (this->bounded_index)
                    return {
                        this->diameter.size(),
                        this->material.size(),
                        this->n_medium.size()
                    };

                if (!this->bounded_index)
                    return {
                        this->diameter.size(),
                        this->index.size(),
                        this->n_medium.size()
                    };

            }

            size_t get_array_size()
            {
                std::vector<size_t> full_shape = this->get_array_shape();
                size_t full_size = 1;
                for (auto e: full_shape)
                    full_size *= e;

                return full_size;
            }
    };
}


namespace CYLINDER
{
    class Set
    {
        public:
            DVector diameter;
            std::vector<complex128> index;
            std::vector<std::vector<complex128>> material;
            DVector n_medium;
            bool bounded_index;
            std::vector<State> States;

            Set(){}
            Set(DVector &diameter, std::vector<std::vector<complex128>> &material, DVector &n_medium)
            : diameter(diameter), material(material), n_medium(n_medium){bounded_index=true;}

            Set(DVector &diameter, std::vector<complex128> &index, DVector &n_medium)
            : diameter(diameter), index(index), n_medium(n_medium){bounded_index=false;}

            State operator[](size_t idx){return this->States[idx];}

            std::vector<size_t> get_array_shape()
            {
                if (this->bounded_index)
                    return {
                        this->diameter.size(),
                        this->material.size(),
                        this->n_medium.size()
                    };

                if (!this->bounded_index)
                    return {
                        this->diameter.size(),
                        this->index.size(),
                        this->n_medium.size()
                    };

            }

            size_t get_array_size()
            {
                std::vector<size_t> full_shape = this->get_array_shape();
                size_t full_size = 1;
                for (auto e: full_shape)
                    full_size *= e;

                return full_size;
            }
    };
}


namespace CORESHELL
{
    class Set
    {
        public:
            DVector
                core_diameter,
                shell_width;

            std::vector<complex128>
                core_index,
                shell_index;

            std::vector<CVector>
                core_material,
                shell_material;

            DVector n_medium;
            bool bounded_core, bounded_shell;
            std::vector<State> States;

            Set(){}

            Set(DVector &core_diameter, DVector &shell_width, CVector &core_index, CVector &shell_index, DVector &n_medium)
            : core_diameter(core_diameter), shell_width(shell_width), core_index(core_index), shell_index(shell_index), n_medium(n_medium)
            {bounded_core=false; bounded_shell=false;}

            Set(DVector &core_diameter, DVector &shell_width, CVector &core_index, std::vector<CVector> &shell_material, DVector &n_medium)
            : core_diameter(core_diameter), shell_width(shell_width), core_index(core_index), shell_material(shell_material), n_medium(n_medium)
            {bounded_core=false; bounded_shell=true;}

            Set(DVector &core_diameter, DVector &shell_width, std::vector<CVector> &core_material, CVector &shell_index, DVector &n_medium)
            : core_diameter(core_diameter), shell_width(shell_width), shell_index(shell_index), core_material(core_material), n_medium(n_medium)
            {bounded_core=true; bounded_shell=false;}

            Set(DVector &core_diameter, DVector &shell_width, std::vector<CVector> &core_material, std::vector<CVector> &shell_material, DVector &n_medium)
            : core_diameter(core_diameter), shell_width(shell_width), core_material(core_material), shell_material(shell_material), n_medium(n_medium)
            {bounded_core=true; bounded_shell=true;}

            State operator[](size_t idx){return this->States[idx];}

            std::vector<size_t> get_array_shape()
            {

                if (this->bounded_core && this->bounded_shell)
                    return {
                        this->core_diameter.size(),
                        this->shell_width.size(),
                        this->core_material.size(),
                        this->shell_material.size(),
                        this->n_medium.size()
                    };

                if (this->bounded_core && !this->bounded_shell)
                    return {
                        this->core_diameter.size(),
                        this->shell_width.size(),
                        this->core_material.size(),
                        this->shell_index.size(),
                        this->n_medium.size()
                    };

                if (!this->bounded_core && this->bounded_shell)
                    return {
                        this->core_diameter.size(),
                        this->shell_width.size(),
                        this->core_index.size(),
                        1,
                        this->n_medium.size()
                    };

                if (!this->bounded_core && !this->bounded_shell)
                    return {
                        this->core_diameter.size(),
                        this->shell_width.size(),
                        this->core_index.size(),
                        this->shell_index.size(),
                        this->n_medium.size()
                    };
            }

            size_t get_array_size()
            {
                std::vector<size_t> full_shape = this->get_array_shape();
                size_t full_size = 1;
                for (auto e: full_shape)
                    full_size *= e;

                return full_size;
            }
    };
}





namespace SOURCE
{
    class Set
    {
        public:
            DVector wavelength;
            std::vector<CVector> jones_vector;
            DVector amplitude;
            std::vector<State> States;

            Set(){}
            Set(DVector &wavelength, std::vector<CVector> &jones_vector, DVector &amplitude)
            : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude)
            {}

            State operator[](size_t idx){return this->States[idx];}

            std::vector<size_t> get_array_shape()
            {
                return {
                    this->wavelength.size(),
                    this->jones_vector.size(),
                    this->amplitude.size(),
                    };
            }

            size_t get_array_size()
            {
                std::vector<size_t> full_shape = this->get_array_shape();
                size_t full_size = 1;
                for (auto e: full_shape)
                    full_size *= e;

                return full_size;
            }
    };
}

namespace DETECTOR
{
    class Set
    {
        public:
            std::vector<CVector> scalar_field;
            DVector NA, phi_offset, gamma_offset, polarization_filter;
            bool coherent, point_coupling;
            std::vector<State> States;

            Set(){}
            Set(std::vector<CVector> &scalar_field, DVector &NA, DVector &phi_offset, DVector &gamma_offset, DVector &polarization_filter, bool    &coherent, bool    &point_coupling)
            : scalar_field(scalar_field), NA(NA), phi_offset(phi_offset), gamma_offset(gamma_offset), polarization_filter(polarization_filter), coherent(coherent), point_coupling(point_coupling)
            {}

            State operator[](size_t idx){return this->States[idx];}

                std::vector<size_t> get_array_shape()
                {
                    return {
                        this->scalar_field.size(),
                        this->NA.size(),
                        this->phi_offset.size(),
                        this->gamma_offset.size(),
                        this->polarization_filter.size()
                        };
                }

                size_t get_array_size()
                {
                    std::vector<size_t> full_shape = this->get_array_shape();
                    size_t full_size = 1;
                    for (auto e: full_shape)
                        full_size *= e;

                    return full_size;
                }
    };
}


#endif
