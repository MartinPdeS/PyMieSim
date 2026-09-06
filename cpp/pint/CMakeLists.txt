# cpp/pint/CMakeLists.txt
add_library(pint STATIC pint.cpp)
target_link_libraries(pint PUBLIC pybind11::pybind11)

pybind11_add_module(interface_pint MODULE interface.cpp)
set_target_properties(interface_pint PROPERTIES OUTPUT_NAME "_pint")
target_link_libraries(interface_pint PRIVATE pint pybind11::module)

install(
    TARGETS interface_pint
    LIBRARY DESTINATION "PyMieSim"
    RUNTIME DESTINATION "PyMieSim"
)
