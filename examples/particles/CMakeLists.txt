file (GLOB CORGI_PARTICLES_SOURCES *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module (corgi_pyprtcls ${CORGI_PARTICLES_SOURCES})
target_link_libraries (corgi_pyprtcls PRIVATE corgi corgi_warnings)

set_target_properties (corgi_pyprtcls PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
)
