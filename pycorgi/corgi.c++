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
PYBIND11_MODULE(corgi, m) {

  /*
    m.attr("Nx")     = conf::Nx;
    m.attr("Ny")     = conf::Ny;
    m.attr("NxCell") = conf::NxCell;
    m.attr("NyCell") = conf::NyCell;
    m.attr("xmin")   = conf::xmin;
    m.attr("xmax")   = conf::xmax;
    m.attr("ymin")   = conf::ymin;
    m.attr("ymax")   = conf::ymax;

    m.def("setSize",       &setSize);
    m.def("setGridLims",   &setGridLims);
    // m.def("printConf",     &printConf);
  */


    // py::class_<corgi::Cell> corgiCell(m, "Cell" );
    //py::class_<Base, std::shared_ptr<Base>>(...);
    // py::class_<Derived, Base, std::shared_ptr<Derived>>(...);
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
        .def("index",                                 &corgi::Cell::index)
        .def("neighs",                                &corgi::Cell::neighs)
        .def("nhood",                                 &corgi::Cell::nhood);


    py::class_<corgi::Node> corgiNode(m, "Node" );
    corgiNode
        .def(py::init<size_t, size_t>())
        .def_readwrite("rank",       &corgi::Node::rank)
        .def_readwrite("Nrank",      &corgi::Node::Nrank)
        .def_readwrite("master",     &corgi::Node::master)

        .def("setGridLims",          &corgi::Node::setGridLims)
        .def("getNx",                &corgi::Node::getNx)
        .def("getNy",                &corgi::Node::getNy)
        .def("getXmin",              &corgi::Node::getXmin)
        .def("getXmax",              &corgi::Node::getXmax)
        .def("getYmin",              &corgi::Node::getYmin)
        .def("getYmax",              &corgi::Node::getYmax)

        .def("getMpiGrid",           &corgi::Node::pyGetMpiGrid)
        .def("setMpiGrid",           &corgi::Node::pySetMpiGrid)

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

        .def("cellId",               &corgi::Node::cellId)
        .def("addCell",              &corgi::Node::addCell)
        .def("getCellIds",           &corgi::Node::getCellIds,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("getCellPtr", py::overload_cast<const uint64_t>(&corgi::Node::getCellPtr),             py::return_value_policy::automatic)
        .def("getCellPtr", py::overload_cast<const size_t, const size_t>(&corgi::Node::getCellPtr), py::return_value_policy::automatic)
        //.def("getCell",    py::overload_cast<const uint64_t>(&corgi::Node::getCell),                py::return_value_policy::reference)
        //.def("getCell",    py::overload_cast<const size_t, const size_t>(&corgi::Node::getCell),    py::return_value_policy::reference);
        .def("getCell",    &corgi::Node::getCell, py::return_value_policy::reference);

        // .def("getCells",             &corgi::Node::getCells,
        //         py::arg("criteria") = std::vector<int>(),
        //         py::arg("sorted") = true)
        // .def("getVirtuals",          &corgi::Node::getVirtuals,
        //         py::arg("criteria") = std::vector<int>(),
        //         py::arg("sorted") = true)

        // .def("isLocal",              &corgi::Node::isLocal)
        // .def("analyzeBoundaryCells", &corgi::Node::analyzeBoundaryCells)


        // // communication wrappers
        // .def_readwrite("send_queue",         &corgi::Node::send_queue)
        // .def_readwrite("send_queue_address", &corgi::Node::send_queue_address)
        // .def("setMpiGrid",           &corgi::Node::setMpiGrid)
        // .def("initMpi",              &corgi::Node::initMpi)
        // .def("bcastMpiGrid",         &corgi::Node::bcastMpiGrid)
        // .def("communicateSendCells", &corgi::Node::communicateSendCells)
        // .def("communicateRecvCells", &corgi::Node::communicateRecvCells)
        // .def("finalizeMpi",          &corgi::Node::finalizeMpi);

}




