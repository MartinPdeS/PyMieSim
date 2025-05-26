# cpp/sets/CMakeLists.txt

set(NAME "scatterer")


add_subdirectory(base_scatterer)     # cpp_base_scatterer
add_subdirectory(sphere)             # cpp_sphere         - interface_sphere
add_subdirectory(cylinder)           # cpp_cylinder       - interface_cylinder
add_subdirectory(coreshell)          # cpp_coreshell      - interface_coreshell


# Create a Python module.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE cpp_sphere cpp_cylinder cpp_coreshell)


install(
    TARGETS interface_${NAME}
    LIBRARY DESTINATION ${PYMIESIM_BIN_DIR}
    ARCHIVE DESTINATION ${PYMIESIM_BIN_DIR}
    RUNTIME DESTINATION ${PYMIESIM_BIN_DIR}
)