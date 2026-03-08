# cpp/sets/CMakeLists.txt

set(NAME "sets")

# Create a shared library for functionality.
add_library(
    "${NAME}"
    STATIC
    "scatterer_set/sphere.cpp"
    "scatterer_set/cylinder.cpp"
    "scatterer_set/core_shell.cpp"
    "source_set/source_set.cpp"
    "detector_set/detector_set.cpp"
    "properties_set/properties_set.cpp"
)

target_link_libraries("${NAME}" PRIVATE pybind11::module sphere cylinder coreshell source full_mesh polarization)
