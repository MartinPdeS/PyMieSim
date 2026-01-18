# cpp/sets/CMakeLists.txt

set(NAME "setup_single")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module OpenMP::OpenMP_CXX cylinder sphere coreshell detector source full_mesh)
