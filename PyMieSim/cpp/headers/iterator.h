#include <iostream>
#include <tuple>
#include <vector>
#include <cassert>


template<typename T_0, typename T_1, typename T_2>
class Iterator {
    public:
        std::vector<T_0>* vector_0;
        std::vector<T_1>* vector_1;
        std::vector<T_2>* vector_2;
        size_t current_index;
        std::vector<size_t> shape;

        Iterator() = default;
        // Constructor
        Iterator(const std::vector<T_0>& input_0, const std::vector<T_1>& input_1, const std::vector<T_2>& input_2)
        : vector_0(&input_0), vector_1(&input_1), vector_2(&input_2), current_index(0){
                this->shape = {this->vector_0->size(), this->vector_1->size(), this->vector_2->size()};
            }

        // Reset the iterator to the beginning
        void reset() {
            current_index = 0;
        }

        // Check if next tuple is available
        bool hasNext() const {
            return current_index < vector_0->size() && current_index < vector_1->size() && current_index < vector_2->size();
        }

        // Get the next tuple
        std::tuple<T_0, T_1, T_2> next() {
            assert(hasNext() && "No more elements to iterate.");
            return {vector_0->at(current_index), vector_1->at(current_index), vector_2->at(current_index++)};
        }
};