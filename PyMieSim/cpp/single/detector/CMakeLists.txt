# PyMieSim/cpp/single/detector/CMakeLists.txt
set(NAME "detector")
set(LIB_NAME "${NAME}_lib")

add_library("${LIB_NAME}" STATIC "detector.cpp")
target_link_libraries("${LIB_NAME}" PUBLIC mode_field_lib scatterer_lib mesh_lib optical_interface_lib)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${LIB_NAME}" pint)

install(
    TARGETS "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
)