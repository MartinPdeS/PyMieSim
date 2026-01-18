#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cmath>


enum class AngleUnit {
    Degree,
    Radian
};

inline double degrees_to_radians(double degrees)
{
    return degrees * (3.141592653589793238462643383279502884 / 180.0);
}

inline double radians_to_degrees(double radians)
{
    return radians * (180.0 / 3.141592653589793238462643383279502884);
}

class BasePolarization {
public:
    virtual ~BasePolarization() = default;
};

class JonesVector : public BasePolarization {
public:
    using Complex = std::complex<double>;
    using Row = std::array<Complex, 2>;
    using Element = std::vector<Row>;

    JonesVector() = default;

    // Python equivalent: np.atleast_2d(element).astype(complex)
    // For C plus plus, accept either a single row (2 complex values) or a list of rows.
    explicit JonesVector(const Row& single_row)
    {
        this->element_.push_back(single_row);
    }

    explicit JonesVector(Element rows)
        : element_(std::move(rows))
    {
        if (this->element_.empty()) {
            throw std::invalid_argument("JonesVector element cannot be empty.");
        }
    }

    const Element& elements() const noexcept
    {
        return this->element_;
    }

    std::size_t number_of_rows() const noexcept
    {
        return this->element_.size();
    }

    bool has_angle() const noexcept
    {
        return this->has_angle_;
    }

    // Valid only if has_angle() is true.
    const std::vector<double>& angles() const
    {
        if (!this->has_angle_) {
            throw std::logic_error("This polarization has no angle attribute.");
        }
        return this->angles_;
    }

    AngleUnit angle_unit() const
    {
        if (!this->has_angle_) {
            throw std::logic_error("This polarization has no angle attribute.");
        }
        return this->angle_unit_;
    }

    // Python equivalent of __add__ stacking behavior.
    // If both have angles, angles are concatenated.
    // If only one has angles, result has no angles (Python deletes angle if missing on either side).
    JonesVector operator+(const JonesVector& other) const
    {
        JonesVector output;

        output.element_.reserve(this->element_.size() + other.element_.size());
        output.element_.insert(output.element_.end(), this->element_.begin(), this->element_.end());
        output.element_.insert(output.element_.end(), other.element_.begin(), other.element_.end());

        if (this->has_angle_ && other.has_angle_) {
            if (this->angle_unit_ != other.angle_unit_) {
                throw std::invalid_argument("Cannot add polarizations with different angle units.");
            }

            output.has_angle_ = true;
            output.angle_unit_ = this->angle_unit_;

            output.angles_.reserve(this->angles_.size() + other.angles_.size());
            output.angles_.insert(output.angles_.end(), this->angles_.begin(), this->angles_.end());
            output.angles_.insert(output.angles_.end(), other.angles_.begin(), other.angles_.end());
        }

        return output;
    }

    // String representation similar to Python __str__:
    // If angle exists, show angles, else show element rows.
    std::string to_string() const
    {
        std::ostringstream stream;

        if (this->has_angle_) {
            stream << "[";
            for (std::size_t i = 0; i < this->angles_.size(); ++i) {
                if (i > 0) {
                    stream << ", ";
                }
                stream << this->angles_[i];
            }
            stream << "]";

            stream << " ";
            stream << (this->angle_unit_ == AngleUnit::Degree ? "degree" : "radian");

            return stream.str();
        }

        stream << "[";
        for (std::size_t i = 0; i < this->element_.size(); ++i) {
            if (i > 0) {
                stream << ", ";
            }

            stream << "[";
            stream << this->element_[i][0] << ", " << this->element_[i][1];
            stream << "]";
        }
        stream << "]";

        return stream.str();
    }

protected:
    Element element_;

    bool has_angle_ = false;
    std::vector<double> angles_;
    AngleUnit angle_unit_ = AngleUnit::Degree;

    void set_angle(std::vector<double> angles, AngleUnit unit)
    {
        if (angles.empty()) {
            throw std::invalid_argument("Angle vector cannot be empty.");
        }
        this->has_angle_ = true;
        this->angles_ = std::move(angles);
        this->angle_unit_ = unit;
    }

    void clear_angle() noexcept
    {
        this->has_angle_ = false;
        this->angles_.clear();
        this->angle_unit_ = AngleUnit::Degree;
    }
};

class RightCircular : public JonesVector {
public:
    RightCircular()
        : JonesVector(Row{Complex(1.0, 0.0), Complex(0.0, 1.0)})
    {
    }
};

class LeftCircular : public JonesVector {
public:
    LeftCircular()
        : JonesVector(Row{Complex(1.0, 0.0), Complex(0.0, -1.0)})
    {
    }
};

class Linear : public JonesVector {
public:
    // Mirrors Python Linear.__post_init__:
    // angle stored in original units, but element is built from degrees.
    //
    // In Python:
    //   self.angle = np.atleast_1d(element.magnitude) * element.units
    //   angle = element.to(degree).magnitude
    //   element[:,0] = cos(deg2rad(angle))
    //   element[:,1] = sin(deg2rad(angle))
    //
    // Here:
    //   angles_input are numeric magnitudes in the given unit
    //   element built using degrees
    explicit Linear(double single_angle, AngleUnit unit)
        : Linear(std::vector<double>{single_angle}, unit)
    {
    }

    explicit Linear(std::vector<double> angles_input, AngleUnit unit)
    {
        if (angles_input.empty()) {
            throw std::invalid_argument("Linear polarization requires at least one angle.");
        }

        // Store angles in original unit, like Python stores angle with units
        this->set_angle(angles_input, unit);

        // Convert to degrees for constructing the Jones vectors
        std::vector<double> angles_degrees;
        angles_degrees.reserve(angles_input.size());

        if (unit == AngleUnit::Degree) {
            angles_degrees = angles_input;
        } else {
            for (double rad : angles_input) {
                angles_degrees.push_back(radians_to_degrees(rad));
            }
        }

        this->element_.resize(angles_degrees.size());

        for (std::size_t i = 0; i < angles_degrees.size(); ++i) {
            const double theta_rad = degrees_to_radians(angles_degrees[i]);
            const double c = std::cos(theta_rad);
            const double s = std::sin(theta_rad);

            this->element_[i][0] = Complex(c, 0.0);
            this->element_[i][1] = Complex(s, 0.0);
        }
    }
};

