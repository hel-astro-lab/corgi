#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "example.h"
#include "../cell.h"




PYBIND11_MODULE(example, m) {

  // --------------------------------------------------
  // Loading cell bindings from corgi library
  py::module::import("corgi");
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");


  /*
    //py::class_<corgi::Cell, std::shared_ptr<corgi::Cell> > corgiCell(m, "Cell" );
    //corgiCell
    py::class_<corgi::Cell >(m, "Cell" )
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

        .def("cellId",               &corgi::Node::cellId)
        .def("addCell",              &corgi::Node::addCell)
        .def("getCellIds",             &corgi::Node::getCellIds,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("getCellPtr", py::overload_cast<const uint64_t>(&corgi::Node::getCellPtr),             py::return_value_policy::reference)
        .def("getCellPtr", py::overload_cast<const size_t, const size_t>(&corgi::Node::getCellPtr), py::return_value_policy::reference);
  */



  // Welsh and Pembroke inherit grid infrastructure from base class and then extend it.
  // NOTE: they must be defined with a special shared container (shared_ptr) in order
  // for the python referencing to work correctly.
  //
  // Example structure from pybind:
  // py::class_<example::Welsh>(m, "Welsh", corgiCell)
  // py::class_<Derived, Base, std::shared_ptr<Derived>>(...);
    

  //py::class_<example::Welsh, corgi::Cell, std::shared_ptr<example::Welsh>>(m, "Welsh")
  py::class_<example::Welsh, std::shared_ptr<example::Welsh>>(m, "Welsh", corgiCell )
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("bark",     &example::Welsh::bark);


  //py::class_<example::Pembroke, corgi::Cell, std::shared_ptr<example::Pembroke>>(m, "Pembroke")
  py::class_<example::Pembroke, std::shared_ptr<example::Pembroke>>(m, "Pembroke", corgiCell)
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("howl",     &example::Pembroke::howl)
    .def("bark",     &example::Pembroke::bark);



  // Vallhund is a compound class build using multiple inheritance.
  // It inherits from Swede and Viking, and hence, we need to pybind them too. 
  // NOTE: All of these parent classes must also have a same special container,
  // i.e., shared_ptr (instead of default unique_ptr)
    
  py::class_<example::Swede,
             std::shared_ptr<example::Swede>>(m, "Swede")
    .def(py::init<int>())
    .def_readwrite("number",  &example::Swede::number)
    .def("fika", &example::Swede::fika);
  

  py::class_<example::Viking,
             std::shared_ptr<example::Viking>>(m, "Viking")
    .def(py::init<>())
    .def("prayForOdin", &example::Viking::prayForOdin);

  py::class_<example::Vallhund, 
             example::Swede, 
             example::Viking, 
             std::shared_ptr<example::Vallhund>>(m, "Vallhund", corgiCell) 
    .def(py::init<size_t, size_t, int, size_t, size_t, int>())
    .def("bark",     &example::Vallhund::bark);



  // --------------------------------------------------
  // Grid bindings

  py::class_<example::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>())
    .def("petShop",   &example::Grid::petShop );


}




