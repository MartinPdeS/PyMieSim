# PyMieSim/cpp/single/source/CMakeLists.txt
set(NAME "source")
set(LIB_NAME "${NAME}_lib")

add_library("${LIB_NAME}" STATIC "${NAME}.cpp")
target_link_libraries("${LIB_NAME}" PUBLIC pint)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${LIB_NAME}")

install(
    TARGETS "${LIB_NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
    ARCHIVE DESTINATION "PyMieSim/single"
)