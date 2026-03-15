# PyMieSim/cpp/experiment/setup/CMakeLists.txt
set(NAME "_setup")

add_library("${NAME}" STATIC "setup.cpp")
set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("${NAME}" PUBLIC OpenMP::OpenMP_CXX source_set scatterer_set detector_set)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "experiment"
    RUNTIME DESTINATION "experiment"
    ARCHIVE DESTINATION "experiment"
)