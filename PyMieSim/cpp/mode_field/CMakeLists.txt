# cpp/detector/CMakeLists.txt

set(NAME "mode_field")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PRIVATE pybind11::module bessel_subroutine)

set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "${NAME}")
