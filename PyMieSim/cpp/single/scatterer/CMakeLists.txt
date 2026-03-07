# PyMieSim/cpp/single/scatterer/CMakeLists.txt
set(NAME "scatterer")

add_library(
    "${NAME}"
    STATIC
    "base_scatterer/base_scatterer.cpp"
    "sphere/sphere.cpp"
    "cylinder/cylinder.cpp"
    "coreshell/coreshell.cpp"
)
set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("${NAME}" PUBLIC pybind11::module pint coordinates mode_field full_mesh fibonacci)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
    RUNTIME DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
    ARCHIVE DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/single"
)
