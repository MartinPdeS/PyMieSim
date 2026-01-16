# cpp/coreshell/CMakeLists.txt

set(NAME "coreshell")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module PUBLIC source full_mesh base_scatterer)

set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "${NAME}")
