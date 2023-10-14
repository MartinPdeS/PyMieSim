#pragma once

    #include "definitions.cpp"
    #include <algorithm>

    double NA2Angle(const double &NA)
    {
        if (NA <= 1.0)
            return asin(NA);

        if (NA >= 1.0)
            return asin(NA-1.0) + PI/2.0;

        return 1.0;
    }

    template <typename T>
    T concatenate_vector(const T &vector_0, const T &vector_1)
    {
        T output_vector = vector_0;
        output_vector.insert( output_vector.end(), vector_1.begin(), vector_1.end() );
        return output_vector;
    }

    template <typename T>
    T concatenate_vector(const T &vector_0, const T &vector_1, const T &vector_2)
    {
        T output_vector = vector_0;
        output_vector.insert( output_vector.end(), vector_1.begin(), vector_1.end() );
        output_vector.insert( output_vector.end(), vector_2.begin(), vector_2.end() );
        return output_vector;
    }

    template <typename T>
    T get_vector_sigma(const std::vector<T> &vector)
    {
        T sigma = 1;
        for (auto e: vector)
          sigma *= e;

        return sigma;
    }

    template <class T>
    T Sum(const std::vector<T>& vector)
    {
        const long unsigned int N = vector.size();
        T sum = 0.;
        for (auto v: vector)
            sum += v;

        return sum;
    }

    template <class T>
    T Sum(const std::vector<T>& vector_0, const std::vector<T>& vector_1)
    {
        const size_t N = vector_0.size();
        T sum = 0.;
        for (size_t iter=0; iter<vector_0.size(); iter++)
            sum += vector_0[iter] * vector_1[iter];

        return sum;
    }


    template <class T>
    void Squared(std::vector<T>& vector)
    {
        for (size_t iter=0; iter<vector.size(); iter++)
            vector[iter] = pow(abs(vector[iter]), 2);
    }

    template <class T>
    std::vector<T> Add(std::vector<T>& vector0, std::vector<T>& vector1)
    {
        std::vector<T> output_vector;
        output_vector.reserve( vector0.size() );

        for (size_t iter=0; iter<vector0.size(); iter++)
            output_vector.push_back( vector0[iter] + vector1[iter] );

        return output_vector;
    }



    void
    Unstructured(uint Sampling, complex128 *array0, complex128 *array1, complex128  scalar, complex128 *output)
    {
        for (size_t p=0; p < Sampling; p++ )
        {
            *output   = scalar * array0[p] * array1[p];
            output++;
        }
    }


    CVector
    Unstructured(CVector &array0, CVector &array1, complex128 &scalar)
    {
        CVector output;
        output.reserve(array1.size());

        for (size_t p=0; p < array1.size(); p++ )
            output.push_back( scalar * array0[p] * array1[p] );

        return output;
    }


    CVector
    Structured_( CVector &STerm, CVector &CosSinTerm, complex128 &scalar)
    {
        CVector output;
        output.reserve(STerm.size() * CosSinTerm.size());

        for (auto S : STerm)
            for (auto Trig : CosSinTerm )
                output.push_back( scalar * S * Trig );

        return output;
    }


    void
    Structured(uint ThetaLength, uint PhiLength, complex128 *array0, complex128 *array1, complex128  scalar, complex128 *output)
    {
        for (uint p=0; p < PhiLength; p++ )
            for (uint t=0; t < ThetaLength; t++ )
            {
                *output   = scalar * array0[p] * array1[t];
                output++;
            }
    }
