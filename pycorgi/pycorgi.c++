#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "../common.h"
#include "../corgi.h"

/*
void setSize(size_t Nx, size_t Ny) {
  conf::Nx = Nx;
  conf::Ny = Ny;
};

void setGridLims(double xmin, double xmax, double ymin, double ymax) {
  conf::xmin = xmin;
  conf::xmax = xmax;
  conf::ymin = ymin;
  conf::ymax = ymax;
};
*/

/*
void printConf() {
  fmt::print("Nx:{}\n", conf::Nx);
  fmt::print("Ny:{}\n", conf::Ny);

  fmt::print("x :{}/{}\n", conf::xmin, conf::xmax);
  fmt::print("y :{}/{}\n", conf::ymin, conf::ymax);
};
*/



// --------------------------------------------------
PYBIND11_MODULE(pycorgi, m) {

    using Node = corgi::Node<2>;
    using Tile = corgi::Tile<2>;

    py::class_<corgi::Communication> corgiComm(m, "Communication");
    corgiComm
        .def_readwrite("owner",                       &corgi::Communication::owner                      )
        .def_readwrite("top_virtual_owner",           &corgi::Communication::top_virtual_owner          )
        .def_readwrite("communications",              &corgi::Communication::communications             ) 
        .def_readwrite("number_of_virtual_neighbors", &corgi::Communication::number_of_virtual_neighbors)
        .def_readwrite("local",                       &corgi::Communication::local                      );


    py::class_<Tile, std::shared_ptr<Tile>> corgiTile(m, "Tile" );
    corgiTile
        .def(py::init<>())

        .def_readwrite("cid",                         &Tile::cid)
        .def_readwrite("communication",               &Tile::communication)
        .def_readwrite("mins",                        &Tile::mins)
        .def_readwrite("maxs",                        &Tile::maxs)
        .def_readwrite("index",                       &Tile::index)
        .def("set_tile_mins",                         &Tile::set_tile_mins)
        .def("set_tile_maxs",                         &Tile::set_tile_maxs)
        //.def("neighs",                                &Tile::neighs)
        .def("neighs", [](Tile &t, size_t i, size_t j){ return t.neighs(i,j); })
        .def("nhood",                                 &Tile::nhood);


    py::class_<Node> corgiNode(m, "Node" );
    corgiNode
        .def(py::init<size_t, size_t>())
        .def_readwrite("rank",       &Node::rank)
        .def_readwrite("Nrank",      &Node::Nrank)
        .def_readwrite("master",     &Node::master)

        .def("getNx",   [](Node &n){ return n.getNx(); })
        .def("getNy",   [](Node &n){ return n.getNy(); })
        .def("getXmin", [](Node &n){ return n.getXmin(); })
        .def("getXmax", [](Node &n){ return n.getXmax(); })
        .def("getYmin", [](Node &n){ return n.getYmin(); })
        .def("getYmax", [](Node &n){ return n.getYmax(); })

        //.def("setGridLims",          &Node::setGridLims)
        .def("setGridLims", [](Node &n, double xmin, double xmax, 
                                        double ymin, double ymax)
            { n.setGridLims({{xmin,ymin}}, {{xmax, ymax}}); })

        .def("getMpiGrid", [](Node &n, const size_t i, const size_t j){ 
            const auto val = n.pyGetMpiGrid(i,j); 
            return val;
            })
        .def("setMpiGrid", [](Node &n, size_t i, size_t j, int val){ n.pySetMpiGrid(val, i, j); })

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

        .def("tileId", [](const Node &n, const size_t i, const size_t j){ return n.id(i,j);})
        .def("addTile",              &Node::addTile)
        .def("getTileIds",           &Node::getTileIds,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("getTilePtr", py::overload_cast<const uint64_t>(&Node::getTilePtr))
        .def("getTile", [](Node &n, const size_t i, const size_t j){ return n.getTileInd(i,j); });

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

}




