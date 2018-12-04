cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

project (corgi)

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_subdirectory(pybind11)
add_subdirectory(mpi4cpp)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/mpi4cpp/include )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/mpi4cpp/tools/optional-lite/include )


##################################################
# dependencies

# hpc stuff 
#find_package (MPI)

#find_package (OpenMP)
#if (OpenMP_FOUND)
##set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif (OpenMP_FOUND)



##################################################
# targets build

#-------------------------------------------------- 
# Warning flags based on the compiler
set (WARNING_FLAGS
    $<$<CXX_COMPILER_ID:GNU>:
        -Wall -Wextra -Wpedantic
    >
    $<$<CXX_COMPILER_ID:Clang>:
        -Wall -Wextra -Wpedantic -Wno-missing-braces 
    >
)

#-------------------------------------------------- 

add_subdirectory(pycorgi)

add_subdirectory(tests)

add_subdirectory(examples/game-of-life)
add_subdirectory(examples/particles)


#-------------------------------------------------- 
# unit tests
#
add_custom_target(check-pycorgi ALL
                  ${PYTHON_EXECUTABLE} -m unittest discover -s ../tests/ -v
                  DEPENDS pycorgi pycorgitest
                  VERBATIM
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
                  )


