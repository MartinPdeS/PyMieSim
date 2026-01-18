# cpp/optical_interface/CMakeLists.txt

set(NAME "optical_interface")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module)

