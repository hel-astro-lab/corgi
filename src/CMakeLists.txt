add_library (corgi)

set_target_properties(corgi PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(corgi PROPERTIES LINKER_LANGUAGE CXX)

target_sources (
  corgi
  PUBLIC
  FILE_SET
  all_headers
  TYPE
  HEADERS
  BASE_DIRS
  "."
  FILES
  ./corgi/cellular_automata.h
  ./corgi/common.h
  ./corgi/communication.h
  ./corgi/corgi.h
  ./corgi/fwd_corgi.h
  ./corgi/internals.h
  ./corgi/tags.h
  ./corgi/tile.h
  ./corgi/toolbox
  ./corgi/geometry/distance.h
  ./corgi/geometry/utilities.h
  ./corgi/toolbox/dataContainer.h
  ./corgi/toolbox/frequency.h
  ./corgi/toolbox/sparse_grid.h
  ./corgi/toolbox/unstable_remove.h
)

target_link_libraries (corgi PUBLIC mpi4cpp PRIVATE corgi_warnings)
