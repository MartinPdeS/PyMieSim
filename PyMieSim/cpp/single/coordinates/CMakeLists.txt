# PyMieSim/cpp/single/coordinates/CMakeLists.txt
set(NAME "coordinates")

add_library("${NAME}" STATIC "${NAME}.cpp")
set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
    RUNTIME DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
    ARCHIVE DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
)