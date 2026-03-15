# PyMieSim/cpp/mesh/CMakeLists.txt
set(NAME "mesh")

add_library("${NAME}" STATIC "fibonacci.cpp" "full_mesh.cpp")
# set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("${NAME}" PUBLIC pint coordinates)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim"
    RUNTIME DESTINATION "PyMieSim"
    ARCHIVE DESTINATION "PyMieSim"
)