#include <vector>
#include <complex>
#include <iostream>
#include <cmath>

constexpr double PI = 3.14159265358979323846;

class BaseSource {
public:
    double wavelength = 0.0;
    std::vector<std::complex<double>> jones_vector;
    double amplitude = 0.0;
    double k = 0.0;

    BaseSource() = default;
    BaseSource(double wavelength, std::vector<std::complex<double>> jones_vector, double amplitude)
    : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude) {
        this->k = 2 * PI / this->wavelength;
    }
};

class Planewave: public BaseSource {
public:
    Planewave() = default;
    Planewave(double wavelength, std::vector<std::complex<double>> jones_vector, double amplitude)
    : BaseSource(wavelength, jones_vector, amplitude) {}
};

class PlanewaveIterator {
private:
    const std::vector<double>& wavelengths;
    const std::vector<std::vector<std::complex<double>>>& jones_vectors;
    const std::vector<double>& amplitudes;
    std::vector<size_t> indices;

public:
    PlanewaveIterator(const std::vector<double>& wavelengths, const std::vector<std::vector<std::complex<double>>>& jones_vectors, const std::vector<double>& amplitudes) :
    wavelengths(wavelengths), jones_vectors(jones_vectors), amplitudes(amplitudes){}

    bool hax_next() const {
        return this->indices[0] < wavelengths.size();
    }

    Planewave next()
    {
        if (!this->hax_next()) {
            throw std::out_of_range("No more combinations.");
        }

        // // Create the Planewave from the current indices
        for (auto amplitude: amplitudes)
            for (auto wavelength: wavelengths)
                for (auto jones_vector: jones_vectors)
                    Planewave test = Planewave(wavelength, jones_vector, amplitude);

        return test;
    }
};

int main() {
    std::vector<double> wavelengths = {600, 700};
    std::vector<std::vector<std::complex<double>>> jonesVectors = {
        {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 1.0)},
        {std::complex<double>(0.5, 0.5), std::complex<double>(0.5, 0.5)}
    };
    std::vector<double> amplitudes = {1.0, 0.5};

    PlanewaveIterator iterator(wavelengths, jonesVectors, amplitudes);

    while (iterator.hax_next()) {
        Planewave planewave = iterator.next();
        std::cout << "Planewave created with wavelength: " << planewave.wavelength << ", amplitude: " << planewave.amplitude << std::endl;
    }

    return 0;
}
