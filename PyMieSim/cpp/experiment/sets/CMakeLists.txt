# cpp/sets/CMakeLists.txt

set(NAME "sets")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PRIVATE pybind11::module sphere cylinder coreshell source full_mesh)
