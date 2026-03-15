# PyMieSim/cpp/pint/CMakeLists.txt
set(NAME "pint")

add_library("${NAME}" STATIC "${NAME}.cpp")
set_target_properties("${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("${NAME}" PUBLIC pybind11::pybind11)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC "${NAME}" pybind11::module)

install(
    TARGETS "${NAME}" "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim/binary"
    RUNTIME DESTINATION "PyMieSim/binary"
    ARCHIVE DESTINATION "PyMieSim/binary"
)



# # PyMieSim/cpp/pint/CMakeLists.txt
# set(NAME "pint")
# set(LIB_NAME "${NAME}_lib")

# add_library("${LIB_NAME}" STATIC "${NAME}.cpp")
# target_link_libraries("${LIB_NAME}" PUBLIC pybind11::module)

# pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
# set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
# target_link_libraries("interface_${NAME}" PUBLIC pybind11::module "${LIB_NAME}")

# install(
#     TARGETS "${LIB_NAME}" "interface_${NAME}"
#     LIBRARY DESTINATION "PyMieSim/binary"
#     RUNTIME DESTINATION "PyMieSim/binary"
#     ARCHIVE DESTINATION "PyMieSim/binary"
# )
