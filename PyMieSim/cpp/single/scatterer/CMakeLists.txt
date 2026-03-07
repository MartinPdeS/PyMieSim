# # cpp/single/CMakeLists.txt

# set(NAME "single")

# add_subdirectory(material)                     # material
# # add_subdirectory(polarization)                 # polarization  - optical_interface
# add_subdirectory(scatterer/base_scatterer)     # base_scatterer
# add_subdirectory(scatterer/sphere)             # sphere         - interface_sphere
# add_subdirectory(scatterer/cylinder)           # cylinder       - interface_cylinder
# add_subdirectory(scatterer/coreshell)          # coreshell      - interface_coreshell

# add_subdirectory(setup_single)
# add_subdirectory(source)
# add_subdirectory(optical_interface)
# add_subdirectory(detector)

# # Create a Python module.
# pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
# set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
# target_link_libraries("interface_${NAME}" PUBLIC detector sphere cylinder coreshell base_scatterer full_mesh mode_field fibonacci)


# set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "interface_${NAME}")


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
