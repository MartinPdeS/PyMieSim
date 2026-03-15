# PyMieSim/cpp/single/scatterer/CMakeLists.txt
set(NAME "scatterer")
set(LIB_NAME "${NAME}_lib")

add_library(
    "${LIB_NAME}"
    STATIC
    "base_scatterer/base_scatterer.cpp"
    "sphere/sphere.cpp"
    "cylinder/cylinder.cpp"
    "coreshell/coreshell.cpp"
)

target_link_libraries("${LIB_NAME}" PUBLIC OpenMP::OpenMP_CXX coordinates_lib mode_field_lib mesh_lib)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${LIB_NAME}" pint)

install(
    TARGETS "${LIB_NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
    ARCHIVE DESTINATION "PyMieSim/single"
)
