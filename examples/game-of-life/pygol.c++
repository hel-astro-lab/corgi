#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "gol.h"
#include "../../toolbox/dataContainer.h"


PYBIND11_MODULE(pyca, m) {


  /// Bind 2D Mesh 
  py::class_<gol::Mesh>(m, "Mesh")
    .def(py::init<int, int>())
    .def_readwrite("Nx",  &gol::Mesh::Nx)
    .def_readwrite("Ny",  &gol::Mesh::Ny)
    .def("__getitem__", [](const gol::Mesh &s, py::tuple indx) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();

        if (i < -s.halo) throw py::index_error();
        if (j < -s.halo) throw py::index_error();

        if (i > s.Nx+s.halo) throw py::index_error();
        if (j > s.Ny+s.halo) throw py::index_error();

        return s(i,j);
      })
    .def("__setitem__", [](gol::Mesh &s, py::tuple indx, int val) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();

        if (i < -s.halo) throw py::index_error();
        if (j < -s.halo) throw py::index_error();

        if (i > s.Nx+s.halo) throw py::index_error();
        if (j > s.Ny+s.halo) throw py::index_error();

        s(i,j) = val;
      });


  // --------------------------------------------------
  // Loading cell bindings from corgi library
  //py::object corgi_tile = (py::object) py::module::import("pycorgi").attr("Tile");

  // cell binding
  py::class_<gol::Tile, 
            corgi::Tile<2>, 
            std::shared_ptr<gol::Tile>
            >(m, "Tile")
    .def(py::init<>())
    .def("add_data",          &gol::Tile::add_data)
    .def("get_data",          &gol::Tile::get_data)
    .def("cycle",             &gol::Tile::cycle)
    .def("update_boundaries", &gol::Tile::update_boundaries);



  // --------------------------------------------------
  // Grid bindings
  //py::object corgi_node = (py::object) py::module::import("corgi").attr("Grid");
  //py::class_<gol::Grid>(m, "Grid", corgi_node)
  //  .def(py::init<size_t, size_t>());
  //  // .def("pet_shop",   &gol::Grid::pet_shop);

  py::class_<gol::Solver>(m, "Solver")
    .def(py::init<>())
    .def("solve", &gol::Solver::solve);

}


