add_library (corgi INTERFACE)

target_sources (
  corgi
  INTERFACE
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
  ./corgi/geometry/distance.h
  ./corgi/toolbox/dataContainer.h
  ./corgi/toolbox/frequency.h
  ./corgi/toolbox/sparse_grid.h
  ./corgi/toolbox/unstable_remove.h
)

target_link_libraries (corgi INTERFACE mpi4cpp corgi_warnings)
