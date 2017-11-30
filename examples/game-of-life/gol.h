#pragma once

#include <vector>

#include "../../cell.h"
#include "../../corgi.h"
#include "../../toolbox/dataContainer.h"



namespace gol {


/// Snapshot patch of a CA simulation
class Patch {

  /// patch dimensions
  size_t Nx;
  size_t Ny;

  /// Internal 2D mesh storing the values
  std::vector<int> mesh;

  /// ctor
  Patch(size_t Nx, size_t Ny);

  /// 2D access operator for values
  int operator () (size_t i, size_t j) {
    size_t indx = Nx*j + i;
    return mesh[indx];
  };

};




/// Small local cellular automata patch
class CellularAutomataCell : public corgi::Cell {

  public:
    CellularAutomataCell(size_t i, size_t j, 
             int o, 
             size_t nx, size_t ny
             ) : corgi::Cell(i, j, o, nx, ny) { }

    ~CellularAutomataCell() { };

    // extending the base class

    datarotators::DataContainer<Patch> data;

    /// Add data to the container
    void addData(Patch m) {
      data.push_back(m);
    }

    /// get current patch
    Patch& getData() {
      return *data.get();
    };




};




/// Simulation grid
class Grid : public corgi::Node {

  public:
    Grid(size_t nx, size_t ny) : corgi::Node(nx, ny) { }

    ~Grid() { };

    // std::string petShop();


};



} // end of namespace gol
