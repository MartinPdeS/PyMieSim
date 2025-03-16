#include <vector>

template <typename T, typename... Ts> T concatenate_vector(const T& first_vector, const Ts&... other_vectors);

template <typename T> T get_vector_sigma(const std::vector<T> &vector);

template <class T> void Squared(std::vector<T>& vector);

template <class T> std::vector<T> Add(std::vector<T>& vector0, std::vector<T>& vector1);

template <class T> std::vector <std::vector<T>> matrix_multiply(std::vector<std::vector<T>> &a, std::vector <std::vector<T>> &b);

std::vector<std::vector<double>> get_rotation_matrix(std::vector<double> rotation_axis, double rotation_angle);
// -
