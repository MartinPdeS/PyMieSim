# cpp/single/scatterer/CMakeLists.txt
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

add_library(scatterer_binding_utils STATIC utils.cpp)
target_link_libraries(scatterer_binding_utils PRIVATE pybind11::pybind11)

pybind11_add_module("interface_${NAME}" MODULE
    interface.cpp
    base_scatterer/interface.cpp
    sphere/interface.cpp
    coreshell/interface.cpp
    cylinder/interface.cpp
)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${LIB_NAME}" pint scatterer_binding_utils)

install(
    TARGETS "${LIB_NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
    ARCHIVE DESTINATION "PyMieSim/single"
)
