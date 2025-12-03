# cpp/cylinder/CMakeLists.txt

set(NAME "cylinder")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PRIVATE pybind11::module PUBLIC source full_mesh base_scatterer)

set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "${NAME}")
