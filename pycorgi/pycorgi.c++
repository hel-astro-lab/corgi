#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "../common.h"
#include "../cell.h"
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


    py::class_<corgi::Cell, std::shared_ptr<corgi::Cell>> corgiCell(m, "Cell" );
    corgiCell
        .def(py::init<size_t, size_t, int, size_t, size_t>())
        .def_readwrite("cid",                         &corgi::Cell::cid)
        .def_readwrite("owner",                       &corgi::Cell::owner)
        .def_readwrite("top_virtual_owner",           &corgi::Cell::top_virtual_owner)
        .def_readwrite("number_of_virtual_neighbors", &corgi::Cell::number_of_virtual_neighbors)
        .def_readwrite("communications",              &corgi::Cell::communications)
        .def_readwrite("i",                           &corgi::Cell::my_i)
        .def_readwrite("j",                           &corgi::Cell::my_j)
        .def_readwrite("local",                       &corgi::Cell::local)
        .def_readwrite("mins",                        &corgi::Cell::mins)
        .def_readwrite("maxs",                        &corgi::Cell::maxs)
        .def("set_tile_mins",                         &corgi::Cell::set_tile_mins)
        .def("set_tile_maxs",                         &corgi::Cell::set_tile_maxs)
        .def("index",                                 &corgi::Cell::index)
        .def("neighs",                                &corgi::Cell::neighs)
        .def("nhood",                                 &corgi::Cell::nhood);


    py::class_<Node> corgiNode(m, "Node" );
    corgiNode
        .def(py::init<size_t, size_t>())
        //.def(py::init([](size_t i, size_t j) { return Node(i,j); }))
        .def_readwrite("rank",       &Node::rank)
        .def_readwrite("Nrank",      &Node::Nrank)
        .def_readwrite("master",     &Node::master)

        .def("getNx",   [](Node &n){ return n.getNx(); })
        .def("getNy",   [](Node &n){ return n.getNy(); })
        .def("getXmin", [](Node &n){ return n.getXmin(); })
        .def("getXmax", [](Node &n){ return n.getXmax(); })
        .def("getYmin", [](Node &n){ return n.getYmin(); })
        .def("getYmax", [](Node &n){ return n.getYmax(); })

        //.def("getNx",                &Node::getNx)
        //.def("getNy",                &Node::getNy)
        //.def("getXmin",              &Node::getXmin)
        //.def("getXmax",              &Node::getXmax)
        //.def("getYmin",              &Node::getYmin)
        //.def("getYmax",              &Node::getYmax)

        //.def("setGridLims",          &Node::setGridLims)
        .def("setGridLims", [](Node &n, double xmin, double xmax, 
                                        double ymin, double ymax)
            { n.setGridLims({{xmin,ymin}}, {{xmax, ymax}}); })

        .def("GetMpiGrid", [](Node &n, const size_t i, const size_t j){ 
            const auto val = n.pyGetMpiGrid(i,j); 
            return val;
            })
        .def("SetMpiGrid", [](Node &n, size_t i, size_t j, int val){ n.pySetMpiGrid(val, i, j); })

        //.def("getMpiGrid",           &Node::pyGetMpiGrid)
        //.def("setMpiGrid",           &Node::pySetMpiGrid)

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

        //.def("cellId",               &Node::cellId)
        .def("cellId", [](const Node &n, const size_t i, const size_t j){ return n.id(i,j);})
        .def("addCell",              &Node::addCell)
        .def("getCellIds",           &Node::getCellIds,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("getCellPtr", py::overload_cast<const uint64_t>(&Node::getCellPtr))
        .def("getCellPtr", py::overload_cast<const size_t, const size_t>(&Node::getCellPtr));

        // .def("getCells",             &Node::getCells,
        //         py::arg("criteria") = std::vector<int>(),
        //         py::arg("sorted") = true)
        // .def("getVirtuals",          &Node::getVirtuals,
        //         py::arg("criteria") = std::vector<int>(),
        //         py::arg("sorted") = true)

        // .def("isLocal",              &ode::isLocal)
        // .def("analyzeBoundaryCells", &ode::analyzeBoundaryCells)


        // // communication wrappers
        // .def_readwrite("send_queue",         &ode::send_queue)
        // .def_readwrite("send_queue_address", &Node::send_queue_address)
        // .def("setMpiGrid",           &Node::setMpiGrid)
        // .def("initMpi",              &Node::initMpi)
        // .def("bcastMpiGrid",         &Node::bcastMpiGrid)
        // .def("communicateSendCells", &Node::communicateSendCells)
        // .def("communicateRecvCells", &Node::communicateRecvCells)
        // .def("finalizeMpi",          &Node::finalizeMpi);

}




