#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "example.h"




PYBIND11_MODULE(example, m) {

  // --------------------------------------------------
  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");

  // Welsh and Pembroke inherit grid infrastructure from base class and then extend it.
  // NOTE: they must be defined with a special shared container (shared_ptr) in order
  // for the python referencing to work correctly.
  //
  // Example structure from pybind:
  // py::class_<example::Welsh>(m, "Welsh", corgiCell)
  // py::class_<Derived, Base, std::shared_ptr<Derived>>(...);
    
  py::class_<example::Welsh, corgi::Cell, std::shared_ptr<example::Welsh>>(m, "Welsh")
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("bark",     &example::Welsh::bark);

  py::class_<example::Pembroke, corgi::Cell, std::shared_ptr<example::Pembroke>>(m, "Pembroke")
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
             std::shared_ptr<example::Vallhund>>(m, "Vallhund") 
    .def(py::init<size_t, size_t, int, size_t, size_t, int>())
    .def("bark",     &example::Vallhund::bark);



  // --------------------------------------------------
  // Grid bindings
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<example::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>())
    .def("petShop",   &example::Grid::petShop);


}




