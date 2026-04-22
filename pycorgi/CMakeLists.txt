file (GLOB CORGI_SOURCES *.c++ *.h)

pybind11_add_module (pycorgi ${CORGI_SOURCES})
target_link_libraries (pycorgi PRIVATE corgi corgi_warnings)

# Setup installation if buildling with pip/scikit-build-core.
if (PROJECT_IS_TOP_LEVEL AND DEFINED SKBUILD_PROJECT_NAME)
    install (TARGETS pycorgi DESTINATION .)
endif ()
