# cpp/sets/CMakeLists.txt

set(NAME "base_scatterer")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module cpp_source cpp_full_mesh)

