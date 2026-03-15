# PyMieSim/cpp/single/detector/CMakeLists.txt
set(NAME "detector")

add_library("${NAME}_lib" STATIC "detector.cpp")
target_link_libraries("${NAME}_lib" PUBLIC pint mode_field scatterer mesh optical_interface)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}_lib")

install(
    TARGETS "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
)