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
            std::vector<double>
                diameter;

            std::vector<complex128>
                index;

            std::vector<std::vector<complex128>>
                material;

            std::vector<double>
                n_medium;

            bool
                bounded_index;

            std::vector<State>
                state_list;

            Set(){}
            Set(
                const std::vector<double> &diameter,
                const std::vector<std::vector<complex128>> &material,
                const std::vector<double> &n_medium)
            :
                diameter(diameter),
                material(material),
                n_medium(n_medium)
            {
                bounded_index = true;
            }

            Set(
                const std::vector<double> &diameter,
                const std::vector<complex128> &index,
                const std::vector<double> &n_medium)
            :
                diameter(diameter),
                index(index),
                n_medium(n_medium)
            {
                bounded_index = false;
            }

            State operator[](const size_t &idx){return this->state_list[idx];}

            std::vector<size_t> get_array_shape() const
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

            size_t get_array_size() const
            {
                std::vector<size_t>
                    full_shape = this->get_array_shape();

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
            std::vector<double>
                diameter,
                n_medium;

            std::vector<complex128>
                index;

            std::vector<std::vector<complex128>>
                material;

            bool
                bounded_index;

            std::vector<State>
                state_list;

            Set(){}
            Set(
                const std::vector<double> &diameter,
                const std::vector<std::vector<complex128>> &material,
                const std::vector<double> &n_medium)
            :
                diameter(diameter),
                material(material),
                n_medium(n_medium)
            {
                bounded_index = true;
            }

            Set(
                const std::vector<double> &diameter,
                const std::vector<complex128> &index,
                const std::vector<double> &n_medium)
            :
                diameter(diameter),
                index(index),
                n_medium(n_medium)
            {
                bounded_index = false;
            }

            State operator[](const size_t &idx){return this->state_list[idx];}

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
            std::vector<double>
                core_diameter,
                shell_width,
                n_medium;

            std::vector<complex128>
                core_index,
                shell_index;

            std::vector<std::vector<complex128>>
                core_material,
                shell_material;

            bool
                bounded_core,
                bounded_shell;

            std::vector<State>
                state_list;

            Set(){}

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<complex128> &core_index,
                const std::vector<complex128> &shell_index,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter),
                shell_width(shell_width),
                core_index(core_index),
                shell_index(shell_index),
                n_medium(n_medium)
            {
                bounded_core = false;
                bounded_shell = false;
            }

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<complex128> &core_index,
                const std::vector<std::vector<complex128>> &shell_material,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter),
                shell_width(shell_width),
                core_index(core_index),
                shell_material(shell_material),
                n_medium(n_medium)
            {
                bounded_core = false;
                bounded_shell = true;
            }

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<std::vector<complex128>> &core_material,
                const std::vector<complex128> &shell_index,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter),
                shell_width(shell_width),
                shell_index(shell_index),
                core_material(core_material),
                n_medium(n_medium)
            {
                bounded_core = true;
                bounded_shell = false;
            }

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<std::vector<complex128>> &core_material,
                const std::vector<std::vector<complex128>> &shell_material,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter),
                shell_width(shell_width),
                core_material(core_material),
                shell_material(shell_material),
                n_medium(n_medium)
            {
                bounded_core = true;
                bounded_shell = true;
            }

            State operator[](const size_t &idx){return this->state_list[idx];}

            std::vector<size_t> get_array_shape() const
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

            size_t get_array_size() const
            {
                std::vector<size_t>
                    full_shape = this->get_array_shape();

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
            std::vector<double>
                wavelength,
                amplitude;

            std::vector<std::vector<complex128>>
                jones_vector;

            std::vector<State>
                state_list;

            Set(){}

            Set(
                const std::vector<double> &wavelength,
                const std::vector<std::vector<complex128>> &jones_vector,
                const std::vector<double> &amplitude)
            :
                wavelength(wavelength),
                jones_vector(jones_vector),
                amplitude(amplitude)
            {}

            State operator[](const size_t &idx){return this->state_list[idx];}

            std::vector<size_t> get_array_shape() const
            {
                return {
                    this->wavelength.size(),
                    this->jones_vector.size(),
                    };
            }

            size_t get_array_size() const
            {
                std::vector<size_t>
                    full_shape = this->get_array_shape();

                size_t
                    full_size = 1;

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
            std::vector<std::vector<complex128>>
                scalar_field;

            std::vector<double>
                NA,
                phi_offset,
                gamma_offset,
                polarization_filter,
                rotation_angle;

            bool
                coherent,
                point_coupling;

            std::vector<State>
                state_list;

            Set(){}
            Set(
                const std::vector<CVector> &scalar_field,
                const std::vector<double> &NA,
                const std::vector<double> &phi_offset,
                const std::vector<double> &gamma_offset,
                const std::vector<double> &polarization_filter,
                const std::vector<double> &rotation_angle,
                const bool &coherent,
                const bool &point_coupling)
            :
                scalar_field(scalar_field),
                NA(NA),
                phi_offset(phi_offset),
                gamma_offset(gamma_offset),
                polarization_filter(polarization_filter),
                rotation_angle(rotation_angle),
                coherent(coherent),
                point_coupling(point_coupling)
            {}

            State operator[](const size_t &idx){return this->state_list[idx];}

            std::vector<size_t> get_array_shape() const
            {
                return {
                    this->scalar_field.size(),
                    this->NA.size(),
                    this->phi_offset.size(),
                    this->gamma_offset.size(),
                    this->polarization_filter.size()
                    };
            }

            size_t get_array_size() const
            {
                std::vector<size_t>
                    full_shape = this->get_array_shape();

                size_t full_size = 1;
                for (auto e: full_shape)
                    full_size *= e;

                return full_size;
            }
    };
}


#endif
