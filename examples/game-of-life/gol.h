#pragma once

#include <vector>

#include "../../cell.h"
#include "../../corgi.h"
#include "../../toolbox/dataContainer.h"



namespace gol {


/// Snapshot patch of a CA simulation
class patch {

  /// patch dimensions
  size_t Nx;
  size_t Ny;

  std::vector<int> mesh;

  patch(size_t Nx, size_t Ny);


};



/// Small general local cellular automata patch
class CellularAutomataCell : public corgi::Cell {

  public:
    CellularAutomataCell(size_t i, size_t j, 
             int o, 
             size_t nx, size_t ny
             ) : corgi::Cell(i, j, o, nx, ny) { }

    ~CellularAutomataCell() { };

    // extend the base class
    // std::string bark();

    datarotators::DataContainer<patch> data;


};




/// Simulation grid
class Grid : public corgi::Node {

  public:
    Grid(size_t nx, size_t ny) : corgi::Node(nx, ny) { }

    ~Grid() { };

    // std::string petShop();


};



} // end of namespace gol
