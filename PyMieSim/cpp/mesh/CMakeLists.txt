# PyMieSim/cpp/mesh/CMakeLists.txt
set(NAME "mesh")
set(LIB_NAME "${NAME}_lib")

add_library("${LIB_NAME}" STATIC fibonacci.cpp full_mesh.cpp)
target_link_libraries("${LIB_NAME}" PUBLIC pint coordinates_lib)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${LIB_NAME}")

install(
    TARGETS "${LIB_NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim"
    RUNTIME DESTINATION "PyMieSim"
    ARCHIVE DESTINATION "PyMieSim"
)