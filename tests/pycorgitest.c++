#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "corgitest.h"


PYBIND11_MODULE(pycorgitest, m) {

  // --------------------------------------------------
  // Loading tile bindings from corgi library
  py::object corgi_tile = (py::object) py::module::import("pycorgi.twoD").attr("Tile");

  // Welsh and Pembroke inherit grid infrastructure from base class and then extend it.
  // NOTE: they must be defined with a special shared container (shared_ptr) in order
  // for the python referencing to work correctly.
  //
  // Example structure from pybind:
  // py::class_<corgitest::Welsh>(m, "Welsh", corgi_tile)
  // py::class_<Derived, Base, std::shared_ptr<Derived>>(...);
    
  py::class_<corgitest::Welsh, corgi::Tile<2>, std::shared_ptr<corgitest::Welsh>>(m, "Welsh")
    .def(py::init<>())
    .def("bark",     &corgitest::Welsh::bark);

  py::class_<corgitest::Pembroke, corgi::Tile<2>, std::shared_ptr<corgitest::Pembroke>>(m, "Pembroke")
    .def(py::init<>())
    .def("howl",     &corgitest::Pembroke::howl)
    .def("bark",     &corgitest::Pembroke::bark);




  // Vallhund is a compound class build using multiple inheritance.
  // It inherits from Swede and Viking, and hence, we need to pybind them too. 
  // NOTE: All of these parent classes must also have a same special container,
  // i.e., shared_ptr (instead of default unique_ptr)
    
  py::class_<corgitest::Swede,
             std::shared_ptr<corgitest::Swede>>(m, "Swede")
    .def(py::init<int>())
    .def_readwrite("number",  &corgitest::Swede::number)
    .def("fika", &corgitest::Swede::fika);
  

  py::class_<corgitest::Viking,
             std::shared_ptr<corgitest::Viking>>(m, "Viking")
    .def(py::init<>())
    .def("pray_for_odin", &corgitest::Viking::pray_for_odin);

  py::class_<corgitest::Vallhund, 
             corgitest::Swede, 
             corgitest::Viking, 
             std::shared_ptr<corgitest::Vallhund>>(m, "Vallhund") 
    .def(py::init<int>())
    .def("bark",     &corgitest::Vallhund::bark);



  // --------------------------------------------------
  // Grid bindings
  //py::object corgi_node = (py::object) py::module::import("pycorgi.twoD").attr("Node");
  //py::class_<corgitest::Grid>(m, "Grid", corgi_node)
  //  .def(py::init<size_t, size_t>())
  //  .def("pet_shop",   &corgitest::Grid::pet_shop);


}




