cmake_minimum_required (VERSION 3.23)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set (CMAKE_CXX_STANDARD 20)
    set (CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

set (CMAKE_CXX_EXTENSIONS OFF)

project (
    corgi
    VERSION 1.0.0
    DESCRIPTION "C++ library for massively parallel computational grids."
    HOMEPAGE_URL "https://github.com/hel-astro-lab/corgi"
    LANGUAGES CXX
)

add_subdirectory (mpi4cpp)

set (PYBIND11_FINDPYTHON ON)
find_package (pybind11 CONFIG REQUIRED)

# Use phony target for handling targets.
add_library (corgi_warnings INTERFACE)
target_compile_options (corgi_warnings
    INTERFACE
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic -Wno-missing-braces>
)

add_subdirectory (src)
add_subdirectory (pycorgi)

if (PROJECT_IS_TOP_LEVEL AND NOT DEFINED SKBUILD_PROJECT_NAME)
    message (STATUS "Enabling tests and examples...")
    include (CTest)
    add_subdirectory (tests)

    message(STATUS "Corgi tests ran with python: ${PYTHON_EXECUTABLE}")

    add_test (NAME check-pycorgi COMMAND
        ${Python_EXECUTABLE} -m unittest discover -s ${PROJECT_SOURCE_DIR}/tests/ -v
    )
    set_property (TEST check-pycorgi
        PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/pycorgi:${CMAKE_BINARY_DIR}/tests
    )

    add_subdirectory(examples/game-of-life)
    add_subdirectory(examples/particles)
endif ()
