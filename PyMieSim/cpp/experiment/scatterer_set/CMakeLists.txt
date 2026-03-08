# PyMieSim/cpp/experiment/scatterer_set/CMakeLists.txt
set(NAME "scatterer_set")

add_library("${NAME}" STATIC "sphere_set.cpp" "cylinder_set.cpp" "core_shell_set.cpp")
set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("${NAME}" PUBLIC scatterer)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${NAME}")

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/experiment"
    RUNTIME DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/experiment"
    ARCHIVE DESTINATION "${CMAKE_SOURCE_DIR}/PyMieSim/experiment"
)