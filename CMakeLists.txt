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
    HOMEPAGE_URL "https://github.com/hel-astro-lab"
    LANGUAGES CXX
)

find_package (Python COMPONENTS Interpreter Development REQUIRED)
add_subdirectory (pybind11)
add_subdirectory (mpi4cpp)

# Use phony target for handling targets.
add_library (corgi_warnings INTERFACE)
target_compile_options (corgi_warnings
    INTERFACE
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic -Wno-missing-braces>
)

add_subdirectory (src)
add_subdirectory (pycorgi)

if (PROJECT_IS_TOP_LEVEL)
    message (STATUS "Enabling tests and examples...")
    add_subdirectory (tests)

    message(STATUS "Corgi tests ran with python: ${PYTHON_EXECUTABLE}")

    add_custom_target(check-pycorgi ALL
        ${Python_EXECUTABLE} -m unittest discover -s ../tests/ -v
        DEPENDS pycorgi pycorgitest
        VERBATIM
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
    )

    add_subdirectory(examples/game-of-life)
    add_subdirectory(examples/particles)
endif ()


