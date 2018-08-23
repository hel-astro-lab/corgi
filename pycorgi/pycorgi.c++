#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "../common.h"
#include "../corgi.h"


template<size_t D>
auto declare_tile(
    py::module &m, 
    const std::string& pyclass_name) 
{

    py::class_<corgi::Tile<D>, std::shared_ptr<corgi::Tile<D> >> corgiTile(m, pyclass_name.c_str());

    corgiTile
        .def(py::init<>())

        .def_readwrite("cid",                         &corgi::Tile<D>::cid)
        .def_readwrite("communication",               &corgi::Tile<D>::communication)
        .def_readwrite("mins",                        &corgi::Tile<D>::mins)
        .def_readwrite("maxs",                        &corgi::Tile<D>::maxs)
        .def_readwrite("index",                       &corgi::Tile<D>::index)
        .def("set_tile_mins",                         &corgi::Tile<D>::set_tile_mins)
        .def("set_tile_maxs",                         &corgi::Tile<D>::set_tile_maxs);
        //.def("nhood",                                 &corgi::Tile<D>::nhood)

    return corgiTile;
}


template<size_t D>
auto declare_node(
    py::module &m, 
    const std::string& pyclass_name) 
{

    py::class_<corgi::Node<D> > corgiNode(m, pyclass_name.c_str());

    corgiNode
        .def_readwrite("rank",       &corgi::Node<D>::rank)
        .def_readwrite("Nrank",      &corgi::Node<D>::Nrank)
        .def_readwrite("master",     &corgi::Node<D>::master)


        /*
        .def("getMpiGrid",              [](SparseGrid<int> &s, const size_t i, const size_t j) {
          // size_t i = indx[0].cast<size_t>();
          // size_t j = indx[1].cast<size_t>();
          return s(i,j);
          })
        .def("setMpiGrid",              [](SparseGrid<int> &s, const size_t i, const size_t j, int val) {
          // size_t i = indx[0].cast<size_t>();
          // size_t j = indx[1].cast<size_t>();
          s(i,j) = val;
          })
         */

        .def("addTile",              &corgi::Node<D>::addTile)
        .def("getTileIds",           &corgi::Node<D>::getTileIds,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("getTilePtr", py::overload_cast<const uint64_t>(&corgi::Node<D>::getTilePtr));

        // .def("getTiles",             &Node::getTiles,
        //         py::arg("criteria") = std::vector<int>(),
        //         py::arg("sorted") = true)
        // .def("getVirtuals",          &Node::getVirtuals,
        //         py::arg("criteria") = std::vector<int>(),
        //         py::arg("sorted") = true)

        // .def("isLocal",              &ode::isLocal)
        // .def("analyzeBoundaryTiles", &ode::analyzeBoundaryTiles)


        // // communication wrappers
        // .def_readwrite("send_queue",         &ode::send_queue)
        // .def_readwrite("send_queue_address", &Node::send_queue_address)
        // .def("setMpiGrid",           &Node::setMpiGrid)
        // .def("initMpi",              &Node::initMpi)
        // .def("bcastMpiGrid",         &Node::bcastMpiGrid)
        // .def("communicateSendTiles", &Node::communicateSendTiles)
        // .def("communicateRecvTiles", &Node::communicateRecvTiles)
        // .def("finalizeMpi",          &Node::finalizeMpi);


  return corgiNode;
}




// --------------------------------------------------
PYBIND11_MODULE(pycorgi, m) {


    py::class_<corgi::Communication> corgiComm(m, "Communication");
    corgiComm
        .def_readwrite("owner",                       &corgi::Communication::owner                      )
        .def_readwrite("top_virtual_owner",           &corgi::Communication::top_virtual_owner          )
        .def_readwrite("communications",              &corgi::Communication::communications             ) 
        .def_readwrite("number_of_virtual_neighbors", &corgi::Communication::number_of_virtual_neighbors)
        .def_readwrite("local",                       &corgi::Communication::local                      );
      

    //--------------------------------------------------

    auto n1 = declare_node<1>(m, "Node1D");

    n1.def(py::init<size_t>());
    n1.def("getNx",   [](corgi::Node<1> &n){ return n.getNx(); })
      .def("getXmin", [](corgi::Node<1> &n){ return n.getXmin(); })
      .def("getXmax", [](corgi::Node<1> &n){ return n.getXmax(); })
      .def("getTile", [](corgi::Node<1> &n, const size_t i){ return n.getTileInd(i); })

      .def("setGridLims", [](corgi::Node<1> &n, double xmin, double xmax)
          { n.setGridLims({{xmin}}, {{xmax}}); })

      .def("getMpiGrid", [](corgi::Node<1> &n, const size_t i){ 
          const auto val = n.pyGetMpiGrid(i); 
          return val;
          })
      .def("setMpiGrid", [](corgi::Node<1> &n, size_t i, int val){ n.pySetMpiGrid(val, i); })
      .def("id", [](const corgi::Node<1> &n, const size_t i){ return n.id(i);});
      

    //--------------------------------------------------
    auto n2 = declare_node<2>(m, "Node2D");

    n2.def(py::init<size_t, size_t>());
    n2.def("getNx",   [](corgi::Node<2> &n){ return n.getNx(); })
      .def("getNy",   [](corgi::Node<2> &n){ return n.getNy(); })

      .def("getXmin", [](corgi::Node<2> &n){ return n.getXmin(); })
      .def("getXmax", [](corgi::Node<2> &n){ return n.getXmax(); })
      .def("getYmin", [](corgi::Node<2> &n){ return n.getYmin(); })
      .def("getYmax", [](corgi::Node<2> &n){ return n.getYmax(); })

      .def("getTile", [](corgi::Node<2> &n, const size_t i, const size_t j){ return n.getTileInd(i,j); })

      .def("setGridLims", [](corgi::Node<2> &n, double xmin, double xmax, 
                                      double ymin, double ymax)
          { n.setGridLims({{xmin,ymin}}, {{xmax, ymax}}); })

      .def("getMpiGrid", [](corgi::Node<2> &n, const size_t i, const size_t j){ 
          const auto val = n.pyGetMpiGrid(i,j); 
          return val;
          })
      .def("setMpiGrid", [](corgi::Node<2> &n, size_t i, size_t j, int val){ n.pySetMpiGrid(val, i, j); })
      .def("id", [](const corgi::Node<2> &n, const size_t i, const size_t j){ return n.id(i,j);});



    //--------------------------------------------------

    auto t1 = declare_tile<1>(m, "Tile1D");
    t1.def("neighs", [](corgi::Tile<1> &t, size_t i ){ return t.neighs(i); });



    auto t2 = declare_tile<2>(m, "Tile2D");
    t2.def("neighs", [](corgi::Tile<2> &t, size_t i, size_t j){ return t.neighs(i,j); });

}




