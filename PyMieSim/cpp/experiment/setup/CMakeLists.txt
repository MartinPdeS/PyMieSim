# PyMieSim/cpp/experiment/setup/CMakeLists.txt
set(NAME "_setup")
set(LIB_NAME "${NAME}_lib")

add_library("${LIB_NAME}" STATIC "setup.cpp")
target_link_libraries("${LIB_NAME}" PUBLIC OpenMP::OpenMP_CXX source_set scatterer_set detector_set)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${LIB_NAME}")

install(
    TARGETS "${LIB_NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/experiment"
    RUNTIME DESTINATION "PyMieSim/experiment"
    ARCHIVE DESTINATION "PyMieSim/experiment"
)