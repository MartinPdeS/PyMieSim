# PyMieSim/cpp/single/optical_interface/CMakeLists.txt
set(NAME "optical_interface")

add_library("${NAME}" STATIC "${NAME}.cpp")
set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
    ARCHIVE DESTINATION "PyMieSim/single"
)