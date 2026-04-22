file (GLOB CORGI_SOURCES *.c++ *.h)

pybind11_add_module (pycorgi ${CORGI_SOURCES})
target_link_libraries (pycorgi PRIVATE corgi corgi_warnings)

if (DEFINED SKBUILD_PROJECT_NAME)
  install (TARGETS pycorgi corgi DESTINATION ${SKBUILD_PROJECT_NAME}
    FILE_SET all_headers DESTINATION ${SKBUILD_HEADERS_DIR}
  )
endif ()
