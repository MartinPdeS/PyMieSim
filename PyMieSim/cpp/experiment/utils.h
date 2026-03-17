#pragma once

#include <memory>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Validation
{
    struct NamedSize
    {
        std::string name;
        std::size_t size;
    };

    template<typename T>
    std::size_t get_object_size(const T& object)
    {
        return static_cast<std::size_t>(object.size());
    }

    template<typename T, std::size_t N>
    std::size_t get_object_size(const T (&)[N])
    {
        return N;
    }

    inline void validate_same_non_zero_size(const std::vector<NamedSize>& objects)
    {
        if (objects.empty()) {
            throw std::invalid_argument("No objects were provided for size validation.");
        }

        const std::size_t reference_size = objects.front().size;

        for (const auto& object : objects) {
            if (object.size == 0) {
                throw std::invalid_argument(
                    "Object '" + object.name + "' has size 0."
                );
            }
        }

        for (const auto& object : objects) {
            if (object.size != reference_size) {
                std::ostringstream message;
                message << "Objects do not have matching sizes: ";

                for (std::size_t index = 0; index < objects.size(); ++index) {
                    if (index > 0) {
                        message << ", ";
                    }

                    message
                        << objects[index].name
                        << "="
                        << objects[index].size;
                }

                throw std::invalid_argument(message.str());
            }
        }
    }
}


