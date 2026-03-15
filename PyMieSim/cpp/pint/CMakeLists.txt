# PyMieSim/cpp/pint/CMakeLists.txt
set(NAME "pint")

add_library("${NAME}" STATIC "${NAME}.cpp")
# set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("${NAME}" PUBLIC pybind11::module)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/binary"
    RUNTIME DESTINATION "PyMieSim/binary"
    ARCHIVE DESTINATION "PyMieSim/binary"
)
