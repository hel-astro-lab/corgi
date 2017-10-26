#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "../corgi.h"



// --------------------------------------------------
PYBIND11_MODULE(corgi, m) {

    m.attr("Nx")     = conf::Nx;
    m.attr("Ny")     = conf::Ny;
    m.attr("NxCell") = conf::NxCell;
    m.attr("NyCell") = conf::NyCell;
    m.attr("xmin")   = conf::xmin;
    m.attr("xmax")   = conf::xmax;
    m.attr("ymin")   = conf::ymin;
    m.attr("ymax")   = conf::ymax;

    py::class_<corgi::Cell>(m, "Cell" )
        .def(py::init<size_t, size_t, int >())
        .def_readwrite("cid",                         &corgi::Cell::cid)
        .def_readwrite("owner",                       &corgi::Cell::owner)
        .def_readwrite("top_virtual_owner",           &corgi::Cell::top_virtual_owner)
        .def_readwrite("number_of_virtual_neighbors", &corgi::Cell::number_of_virtual_neighbors)
        .def_readwrite("communications",              &corgi::Cell::communications)
        .def_readwrite("i",                           &corgi::Cell::i)
        .def_readwrite("j",                           &corgi::Cell::j)
        .def_readwrite("local",                       &corgi::Cell::local)
        .def("index",                                 &corgi::Cell::index)
        .def("neighs",                                &corgi::Cell::neighs)
        .def("nhood",                                 &corgi::Cell::nhood);


    py::class_<corgi::Node>(m, "Node" )
        .def(py::init<>())
        .def_readwrite("rank",    &corgi::Node::rank)
        .def_readwrite("Nrank",   &corgi::Node::Nrank)
        .def_readwrite("master",  &corgi::Node::master)
        .def("mpiGrid",           &corgi::Node::mpiGrid)
        .def("is_local",          &corgi::Node::is_local)
        .def("cell_id",           &corgi::Node::cell_id)
        .def("add_local_cell",    &corgi::Node::add_local_cell)
        .def("get_cell",          &corgi::Node::get_cell)
        .def("get_cells",         &corgi::Node::get_cells,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("get_virtuals",      &corgi::Node::get_virtuals,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("analyze_boundary_cells", &corgi::Node::analyze_boundary_cells)

        // communication wrappers
        .def_readwrite("send_queue",         &corgi::Node::send_queue)
        .def_readwrite("send_queue_address", &corgi::Node::send_queue_address)
        .def("set_mpiGrid",                  &corgi::Node::set_mpiGrid)
        .def("init_mpi",                     &corgi::Node::init_mpi)
        .def("bcast_mpiGrid",                &corgi::Node::bcast_mpiGrid)
        .def("communicate_send_cells",       &corgi::Node::communicate_send_cells)
        .def("communicate_recv_cells",       &corgi::Node::communicate_recv_cells)
        .def("finalize_mpi",                 &corgi::Node::finalize_mpi);



}
