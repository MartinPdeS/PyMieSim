#include <pybind11/embed.h>
#include <experiment/setup/setup.h>
#include <utils/numpy_interface.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void test_numpy_copies() {
    std::vector<double> values{1., 2., 3., 4.};
    auto array = vector_to_numpy_copy(values, {2, 2});
    require(array.owndata(), "NumPy copy must own its data");
    require(array.ndim() == 2 && array.shape(1) == 2, "Shape must be preserved");
    values[0] = 99.;
    require(array.at(0, 0) == 1., "Copy must not alias its input");
    array.mutable_at(0, 1) = 88.;
    require(values[1] == 2., "Writing to a copy must not mutate its input");
    values.clear();
    values.shrink_to_fit();
    require(array.at(1, 1) == 4., "Array must outlive input storage");
    auto empty = vector_to_numpy_copy(values);
    require(empty.size() == 0 && empty.ndim() == 1, "Empty vector must remain one-dimensional");
    auto scalar = vector_to_numpy_copy(std::vector<double>{7.}, {});
    require(scalar.ndim() == 0, "Explicit empty shape must produce a scalar");
    for (const auto& shape : {std::vector<size_t>{2}, std::vector<size_t>{std::numeric_limits<size_t>::max(), 2}}) {
        bool rejected = false;
        try { vector_to_numpy_copy(std::vector<double>{1.}, shape); }
        catch (const py::value_error&) { rejected = true; }
        require(rejected, "Invalid or overflowing shapes must be rejected");
    }
}

class FailingSphereSet : public SphereSet {
public:
    explicit FailingSphereSet(const SphereSet& source) : SphereSet(source) {}
    std::shared_ptr<BaseScatterer> get_scatterer_by_index(size_t index) const override {
        if (index == 1) throw std::invalid_argument("injected scatterer initialization failure");
        return SphereSet::get_scatterer_by_index(index);
    }
};

void test_farfield_errors() {
    auto sources = std::make_shared<PlaneWaveSourceSet>(
        std::vector<double>{600e-9, 700e-9}, PolarizationSet(0.), std::vector<double>{1.}, false);
    auto scatterers = std::make_shared<SphereSet>(
        std::vector<double>{100e-9, 200e-9}, MaterialSet(std::vector<complex128>{1.4}),
        MediumSet(std::vector<double>{1.}), false);
    Setup setup(scatterers, sources, nullptr);
    FibonacciMesh mesh(16, 1., 0., 0., 0., 0.);
    const auto [expected, shape] = setup.get_farfields(mesh, 1.);
    require(expected.size() == 4 * 2 * 16, "Far-field dimensions are incorrect");
    require(std::any_of(expected.begin(), expected.end(), [](auto x) { return std::abs(x) > 0.; }),
            "Valid fields must not be all zero");
    // A malformed native mesh previously returned a valid-looking zero array.
    ++mesh.sampling;
    bool rejected = false;
    try { setup.get_farfields(mesh, 1.); }
    catch (const std::runtime_error& error) {
        const std::string message = error.what();
        rejected = message.find("configuration 0") != std::string::npos &&
                   message.find("expected 17") != std::string::npos &&
                   message.find("16 phi and 16 theta") != std::string::npos;
    }
    require(rejected, "Malformed fields must raise a descriptive exception after the parallel loop");
    --mesh.sampling;
    require(std::get<0>(setup.get_farfields(mesh, 1.)) == expected,
            "A failed sweep must not corrupt subsequent results");
    Setup failing_setup(std::make_shared<FailingSphereSet>(*scatterers), sources, nullptr);
    rejected = false;
    try { failing_setup.get_farfields(mesh, 1.); }
    catch (const std::invalid_argument& error) {
        rejected = std::string(error.what()) == "injected scatterer initialization failure";
    }
    require(rejected, "Worker exceptions must retain their type and message on the calling thread");
}

int main() {
    py::scoped_interpreter interpreter{};
    try {
        test_numpy_copies();
        test_farfield_errors();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
