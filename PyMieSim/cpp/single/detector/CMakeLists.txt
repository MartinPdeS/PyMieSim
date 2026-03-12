# PyMieSim/cpp/single/detector/CMakeLists.txt
set(NAME "detector")

add_library("${NAME}" STATIC "detector.cpp")
set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("${NAME}" PUBLIC pint mode_field scatterer mesh optical_interface)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
    RUNTIME DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
    ARCHIVE DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
)