# cpp/sphere/CMakeLists.txt

set(NAME "source")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PRIVATE pybind11::module)

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "${NAME}")


install(
    TARGETS ${NAME} interface_${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)