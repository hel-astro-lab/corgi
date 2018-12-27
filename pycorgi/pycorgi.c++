#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "tuple"

#include "../common.h"
#include "../corgi.h"



PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>, true)


template<size_t D>
auto declare_tile(
    py::module &m, 
    const std::string& pyclass_name) 
{

    py::class_<corgi::Tile<D>, 
              std::shared_ptr<corgi::Tile<D> >> corgi_tile(m, pyclass_name.c_str());

    corgi_tile
        .def(py::init<>())

        .def_readwrite("cid",           &corgi::Tile<D>::cid)
        .def_readwrite("communication", &corgi::Tile<D>::communication)
        .def_readwrite("mins",          &corgi::Tile<D>::mins)
        .def_readwrite("maxs",          &corgi::Tile<D>::maxs)
        .def_readwrite("index",         &corgi::Tile<D>::index)
        .def("set_tile_mins",           &corgi::Tile<D>::set_tile_mins)
        .def("set_tile_maxs",           &corgi::Tile<D>::set_tile_maxs)
        .def("load_metainfo",           &corgi::Tile<D>::load_metainfo)
        .def("nhood",                   &corgi::Tile<D>::nhood);

    return corgi_tile;
}


template<size_t D>
auto declare_node(
    py::module &m, 
    const std::string& pyclass_name) 
{

    py::class_<corgi::Node<D> > corgi_node(m, pyclass_name.c_str());

    corgi_node
        .def("rank",      [](corgi::Node<D>& n) { return n.comm.rank(); })
        .def("size",      [](corgi::Node<D>& n) { return n.comm.size(); })
        .def("master",    [](corgi::Node<D>& n) { return n.comm.rank() == 0; })
        .def("add_tile",              &corgi::Node<D>::add_tile, py::keep_alive<1,2>())
        .def("get_tile_ids",           &corgi::Node<D>::get_tile_ids,
                py::arg("sorted") = true)
        .def("get_tile", (std::shared_ptr<corgi::Tile<D>> (corgi::Node<D>::*)(const uint64_t)) &corgi::Node<D>::get_tileptr)

        .def("get_local_tiles",             &corgi::Node<D>::get_local_tiles,
                 py::arg("sorted") = true)
        .def("get_virtual_tiles",          &corgi::Node<D>::get_virtuals,
                 py::arg("sorted") = true)
        .def("get_boundary_tiles",          &corgi::Node<D>::get_boundary_tiles,
                 py::arg("sorted") = true)

        .def("is_local",              &corgi::Node<D>::is_local)
        .def("analyze_boundaries", &corgi::Node<D>::analyze_boundaries)

        // // communication wrappers
        .def("send_tile",               &corgi::Node<D>::send_tile)
        .def("recv_tile",               &corgi::Node<D>::recv_tile)
        .def_readwrite("send_queue",         &corgi::Node<D>::send_queue)
        .def_readwrite("send_queue_address", &corgi::Node<D>::send_queue_address)
        .def("bcast_mpi_grid",          &corgi::Node<D>::bcast_mpi_grid)
        .def("send_tiles",              &corgi::Node<D>::send_tiles)
        .def("recv_tiles",              &corgi::Node<D>::recv_tiles)
        .def("send_data",               &corgi::Node<D>::send_data)
        .def("recv_data",               &corgi::Node<D>::recv_data)
        .def("wait_data",               &corgi::Node<D>::wait_data)

        // adoption routines
        .def("adopt",                   &corgi::Node<D>::adopt)
        .def("adoption_council",        &corgi::Node<D>::adoption_council);


  return corgi_node;
}




// --------------------------------------------------
PYBIND11_MODULE(pycorgi, m_base) {

    //--------------------------------------------------
    // dimension independent bindings

    py::class_<corgi::Communication> corgi_comm(m_base, "Communication");
    corgi_comm
        .def_readwrite("cid",                         &corgi::Communication::cid                        )
        .def_readwrite("indices",                     &corgi::Communication::indices                    )
        .def_readwrite("owner",                       &corgi::Communication::owner                      )
        .def_readwrite("top_virtual_owner",           &corgi::Communication::top_virtual_owner          )
        .def_readwrite("communications",              &corgi::Communication::communications             ) 
        .def_readwrite("number_of_virtual_neighbors", &corgi::Communication::number_of_virtual_neighbors)
        .def_readwrite("mins",                        &corgi::Communication::mins                       )
        .def_readwrite("maxs",                        &corgi::Communication::maxs                       )
        .def_readwrite("local",                       &corgi::Communication::local                      )
        .def_readwrite("virtual_owners",              &corgi::Communication::virtual_owners                      );
      

    //--------------------------------------------------
    // 1D
      
    py::module m_1d = m_base.def_submodule("oneD", "1D specializations");

    auto n1 = declare_node<1>(m_1d, "Node");

    n1.def(py::init<size_t>());

    // make constructor symmetric and assume 3D shape
    n1.def(py::init( [](size_t nx, size_t /*ny*/) {
      //assert(ny == 1);

      return new corgi::Node<1>(nx);}));
    n1.def(py::init( [](size_t nx, size_t /*ny*/, size_t /*nz*/) {
      // assert(ny == 1);
      // assert(nz == 1);
      return new corgi::Node<1>(nx);}));

    //n1.def(py::init([](size_t i, size_t j) {
    //    return std::unique_ptr<Example>(new Example(arg));
    n1.def("get_tile", [](corgi::Node<1> &n, size_t i, size_t /*j*/){
        return n.get_tileptr_ind(i); });
    n1.def("get_tile", [](corgi::Node<1> &n, size_t i, size_t /*j*/, size_t /*k*/){
        return n.get_tileptr_ind(i); });

    n1.def("get_Nx",   [](corgi::Node<1> &n){ return n.get_Nx(); })
      .def("get_xmin", [](corgi::Node<1> &n){ return n.get_xmin(); })
      .def("get_xmax", [](corgi::Node<1> &n){ return n.get_xmax(); });
      //.def("get_tile", [](corgi::Node<1> &n, size_t i){ 
      //    return n.get_tileptr( std::make_tuple<size_t>(i) ); })

    n1.def("get_Ny",   [](corgi::Node<1> &){ return 1; })
      .def("get_Nz",   [](corgi::Node<1> &){ return 1; });

    n1.def("get_ymin", [](corgi::Node<1> & ){ return 0.0; })
      .def("get_ymax", [](corgi::Node<1> & ){ return 1.0; });

    n1
      .def("set_grid_lims", [](corgi::Node<1> &n, double xmin, double xmax)
          { n.set_grid_lims({{xmin}}, {{xmax}}); });

    n1
      .def("set_grid_lims", [](corgi::Node<1> &n, double xmin, double xmax,
                                                double /*ymin*/, double /*ymax*/
            )
          { n.set_grid_lims({{xmin}}, {{xmax}}); });

    n1
      .def("get_mpi_grid", [](corgi::Node<1> &n, const size_t i){ 
          const auto val = n.py_get_mpi_grid(i); 
          return val;
          })
      .def("get_mpi_grid", [](corgi::Node<1> &n, const size_t i, const size_t ){ 
          const auto val = n.py_get_mpi_grid(i); 
          return val;
          })
      .def("set_mpi_grid", [](corgi::Node<1> &n, size_t i, int val){ n.py_set_mpi_grid(val, i); })
      .def("set_mpi_grid", [](corgi::Node<1> &n, size_t i, size_t /*j*/, int val){ n.py_set_mpi_grid(val, i); })
      .def("id", [](const corgi::Node<1> &n, const size_t i){ return n.id(i);});
      

    //--------------------------------------------------
    auto t1 = declare_tile<1>(m_1d, "Tile");
    t1.def("neighs", [](corgi::Tile<1> &t, size_t i ){ return t.neighs(i); });




    //--------------------------------------------------
    // 2D
    py::module m_2d = m_base.def_submodule("twoD", "2D specializations");


    auto n2 = declare_node<2>(m_2d, "Node");
    n2.def(py::init<size_t, size_t>());

    n2.def(py::init( [](size_t nx, size_t ny, size_t /*nz*/) {
      //assert(nz == 1);
      return new corgi::Node<2>(nx, ny);}));

    n2.def("get_Nx",   [](corgi::Node<2> &n){ return n.get_Nx(); })
      .def("get_Ny",   [](corgi::Node<2> &n){ return n.get_Ny(); });

    n2.def("get_Nz",   [](corgi::Node<2> &){  return 1; });

    n2.def("get_xmin", [](corgi::Node<2> &n){ return n.get_xmin(); })
      .def("get_xmax", [](corgi::Node<2> &n){ return n.get_xmax(); })
      .def("get_ymin", [](corgi::Node<2> &n){ return n.get_ymin(); })
      .def("get_ymax", [](corgi::Node<2> &n){ return n.get_ymax(); });


    n2.def("get_tile", [](corgi::Node<2> &n, size_t i, size_t j){ 
          return n.get_tileptr( std::make_tuple(i,j) ); })
      .def("get_tile", [](corgi::Node<2> &n, size_t i, size_t j, size_t /*k*/){
        return n.get_tileptr_ind(i,j); })
      .def("set_grid_lims", [](corgi::Node<2> &n, 
            double xmin, double xmax, 
            double ymin, double ymax)
          { n.set_grid_lims({{xmin,ymin}}, {{xmax, ymax}}); })

      .def("get_mpi_grid", [](corgi::Node<2> &n, const size_t i, const size_t j){ 
          const auto val = n.py_get_mpi_grid(i,j); 
          return val;
          })
      .def("set_mpi_grid", [](corgi::Node<2> &n, size_t i, size_t j, int val){ n.py_set_mpi_grid(val, i, j); })
      .def("id", [](const corgi::Node<2> &n, const size_t i, const size_t j){ return n.id(i,j);});


    auto t2 = declare_tile<2>(m_2d, "Tile");
    t2.def("neighs", [](corgi::Tile<2> &t, size_t i, size_t j){ return t.neighs(i,j); });


    //--------------------------------------------------
    // TODO: 3D



}




