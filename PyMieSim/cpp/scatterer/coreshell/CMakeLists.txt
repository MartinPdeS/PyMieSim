# cpp/coreshell/CMakeLists.txt

set(NAME "coreshell")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module ZBessel PUBLIC cpp_source cpp_full_mesh cpp_base_scatterer)
