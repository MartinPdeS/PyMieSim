# PyMieSim/cpp/single/setup/CMakeLists.txt
set(NAME "setup_single")
add_library("${NAME}" STATIC "setup.cpp")

target_link_libraries("${NAME}" PUBLIC scatterer_lib source_lib detector_lib)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME setup)
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}" pint)

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
    ARCHIVE DESTINATION "PyMieSim/single"
)